%
% function [xI50, Hill, Einf, Area, r2, EC50, fit_final, p, log, xI50_interval] = ...
%     ICcurve_fit(Conc, Growth, fit_type, opt)
%
%   Inputs:
%   ---------
%
%   Conc:   expected in uM, NOT logged
%   Growth: normalized to control (=1) fit_type
%   fit_type:   'GI50' or 'IC50' (default)
%
%   options:
%       - plotting: plot the curves [false]
%       - priors:   seed for the fitting
%       - ranges:   range for the fitting
%       - fitting:  average/individual for replicates [average]
%       - pcutoff:  cutoff for F-test [0.05]
%       - capped:   cap points with enhanced growth [true]
%       - extrapolrange: maximum fold to extrapolate xI50 [10]
%
%   Outputs:
%   ----------
%   concentrations are output in uM
%
%   IC50 => Growth is relative to control (end/ctrl)
%   GI50 => Growth is relative to growth of the control:
%               (end-day0)/(ctrl-day0)
%           Growth should not comprise the control value; each replicate is
%           a column
%
%   Area => sum of (1-cell count) devided by range --> average per order of
%   magnitude.
%

function [xI50, Hill, Einf, Area, r2, EC50, fit_final, p, log] = ...
    ICcurve_fit(Conc, Growth, fit_type, opt)


% parameters : E0 Emax EC50 HS (all in uM)
priors = [1 .1 median(Conc) 2];

ranges = [
    .975 1.025  %E0
    0 1    %Emax
    max(min(Conc)*1e-4,1e-7) min(max(Conc)*1e2, 1e3)  %E50
    .1 5    % HS
    ]';

plotting = 0;
fitting = 'average';
pcutoff = .05;
capped = true;
extrapolrange = 10;

if ~exist('fit_type','var') || isempty(fit_type)
    fit_type = 'IC50';
end

if exist('opt','var')
    fields = {'plotting', 'priors', 'ranges', 'fitting', 'pcutoff' 'capped' 'extrapolrange'};
    for field = fields
        if isfield(opt,field{:})
            eval([field{:} ' = opt.' field{:} ';'])
        end
    end
end

Conc = ToRow(Conc);
if isvector(Growth)
    Growth=ToColumn(Growth);
    fitting = 'average';
end
assert(all([1 size(Growth,1)]==size(Conc)))


switch fitting
    case 'average'
        g = mean(Growth,2)';
    case 'individual'
        for i=1:size(Growth,2)
            [xI50(i), Hill(i), Einf(i), Area(i), r2(i), EC50(i), fit_final{i}, p(i), log{i}] = ...
                ICcurve_fit(Conc, Growth(:,i)', fit_type, opt);
        end
        return
        
    otherwise
        error('wrong input for ''fitting''')
        
end

% 'IC50'
switch fit_type
    case 'IC50'
        ranges(1,2) = 0; % lowest Emax
    case 'GI50'
        ranges(1,2) = -Inf; % lowest Emax
end


% remove the case of enhanced proliferation to avoid failure of F-test
if capped 
    Conc = Conc(~isnan(g));
    g = g(~isnan(g));
    g = min(g, ranges(2,1));
end

Npara = 4; % N of parameters in the growth curve
[fit_res, gof] = sigmoidal_fit(Conc,g);
[fit_res_flat, gof_flat] = flat_fit(Conc,g);

% F-test for the models
RSS2 = gof.sse;
RSS1 = gof_flat.sse;

df1 = (Npara -1);
df2 = (length(g) -Npara);
F = ( (RSS1-RSS2)/df1 )/( RSS2/df2 );
p = 1-fcdf(F, df1, df2);

xc = 10.^(log10(min(Conc))-1:.05:log10(max(Conc))+1);
r2 = gof.rsquare;

if p>=pcutoff
    xI50 = +Inf;
    EC50 = +Inf;
    Hill = 0;
    Einf = 1- mean(g(end-[2 1 0]));
    log = ['**** USING LINEAR FIT **** r2= ' num2str(gof_flat.rsquare,'%.2f')];
    fit_final = fit_res_flat;
    
    Area = length(g)*(1-fit_res_flat.b);
    
else
    log = ['r2 = ' num2str(gof.rsquare,'%.2f')];
    
    fit_growth = fit_res(xc);
    
    Einf = 1- fit_res.b;
    Hill = fit_res.d;
    fit_final = fit_res;
    
    
    EC50 = fit_res.c;
    Area = sum(1-g);
    
    xI50 = fit_res.c*((((fit_res.a-fit_res.b)/(.5-fit_res.b))-1)^(1/fit_res.d));
    if any(fit_growth<.5) && any(fit_growth>.5) % inter/extrapolation is fine
        log = [log '\t' fit_type ' in the range of data']; 
    elseif all(fit_growth>.5)
        if xI50>extrapolrange*max(Conc) || imag(xI50)~=0
            xI50 = Inf;
            log = [log '\t' fit_type '>' num2str(extrapolrange) '*max(Conc) --> +Inf'];
        end
    elseif all(fit_growth<.5) || imag(xI50)~=0
        if xI50<min(Conc)/extrapolrange
            xI50 = -Inf;
            log = [log '\t' fit_type '<min(Conc)/' num2str(extrapolrange) ' --> -Inf'];
        end
    else
        xI50 = NaN;
        warning(['undefined ' fit_type])
        log = [log '\tundefined ' fit_type ' --> NaN'];
    end
    
end


if plotting
    
    errorbar(log10(Conc), mean(Growth,2), std(Growth,[],2) ,'.k-');
    
    hold on
    plot(log10(xc), fit_res(xc),'-r')
    
    plot(log10(EC50)*[1 1], [0 fit_res.b+(fit_res.a-fit_res.b)/2], '.b-')
    plot(log10(Conc([1 end]))+[1 0], [1 1]*fit_res.b, '.b-')
    
    plot(log10(xI50)*[1 1], [0 .5], '.b-')
    
    
    if p>=.05
        plot(log10(xc), fit_res_flat(xc),'-g');
    end
    
    title(sprintf('r^2 = %.3f', r2))
    xlim(log10([min(Conc)/extrapolrange extrapolrange*max(Conc)]))
    ylim([min(min(Growth(:)), max(-.5,1-Einf)) max(Growth(:))*1.05])
end


    function [fit_result, gof2] = sigmoidal_fit(doses, response)
        fitopt = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',ranges(1,:),...
            'Upper',ranges(2,:),...
            'Startpoint',priors);
        f = fittype('b + (a-b) ./ ( 1 + (x/c).^d)','options',fitopt);
        [fit_result,gof2] = fit(doses', response',f);
    end

    function [fit_result, gof2] = flat_fit(doses, response)
        fitopt = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',ranges(1,2),...   % min Emax
            'Upper',ranges(2,1),...   % max E0
            'Startpoint',priors(1));
        f = fittype('b+0*x','options',fitopt);
        [fit_result,gof2] = fit(doses', response',f,'Robust','on');
    end

end