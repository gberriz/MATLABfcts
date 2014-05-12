%
% function [xI50, Hill, Emax, Area, r2, fit_final, p, log] = ...
%     ICcurve_fit(Conc, Growth, fit_type, opt)
%
%   IC50 => Growth is relative to control (end/ctrl)
%   GI50 => Growth is relative to growth of the control:
%               (end-day0)/(ctrl-day0)
%           Growth should not comprise the control value; each replicate is
%           a column
%
%   Area => sum of (1-cell count)
%
%

function [xI50, Hill, Emax, Area, r2, fit_final, p, log] = ...
    ICcurve_fit(Conc, Growth, fit_type, opt)


% parameters : E0 Emax EC50 HS (all in uM, in log10 domain)
priors = [1 .1 1 2];

ranges = [
    .975 1.025  %E0
    0 1    %Emax
    10.^[-4.5 1.5]  %E50
    .1 4    % HS
    ]';

plotting = 0;
fitting = 'average';
pcutoff = .05;

if ~exist('fit_type','var') || isempty(fit_type)
    fit_type = 'IC50';
end

if exist('opt','var')
    fields = {'plotting', 'priors', 'ranges', 'fitting', 'pcutoff'};
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
            [xI50(i), Hill(i), Emax(i), Area(i), r2(i), fit_final{i}, p(i), log{i}] = ...
                ICcurve_fit(Conc, Growth(:,i)', fit_type, opt);
        end
        return
end

% 'IC50'
switch fit_type
    case 'IC50'
        ranges(1,2) = 0; % lowest Emax
    case 'GI50'
        ranges(1,2) = -Inf; % lowest Emax
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
    Hill = 0;
    Emax = fit_res_flat.b;
    log = ['**** USING LINEAR FIT **** r2= ' num2str(gof_flat.rsquare,'%.2f')];
    fit_final = fit_res_flat;
    
    Area = length(g)*(1-fit_res_flat.b);
    
else
    log = ['r2 = ' num2str(gof.rsquare,'%.2f')];
    
    fit_growth = fit_res(xc);
    
    Emax = fit_res.b;
    Hill = fit_res.d;
    fit_final = fit_res;
    
    
    Area = sum(1-g);
    
    switch fit_type
        
        case {'IC50' 'x50'}
            xI50 = fit_res.c;
        case 'GI50'
            if any(fit_growth<.5) && any(fit_growth>.5) % interpolation
                xI50 = fit_res.c*((((fit_res.a-fit_res.b)/(.5-fit_res.b))-1)^(1/fit_res.d));
            elseif all(fit_growth>.5)
                xI50 = Inf;
                log = [log '\tGI50>10*max(Conc) --> +Inf'];
            elseif all(fit_growth<.5)
                xI50 = -Inf;
                log = [log '\tGI50<min(Conc)/10 --> -Inf'];
            else
                xI50 = NaN;
                warning('undefined GI50')
                log = [log '\tundefined GI50 --> NaN'];
            end
            
    end
end


if plotting
    
    if exist('RelativeGrowth','var')
        errorbar(log10(Conc), mean(RelativeGrowth,1), std(RelativeGrowth,[],1) ,'.k-');
    end
    hold on
    plot(log10(xc), fit_res(xc),'-r')
    
    if ismember(fit_type,[0 2])
        plot(log10(IC50)*[1 1], [0 Emax+(fit_res.a-Emax)/2], '.b-')
        plot(log10(Conc([1 end]))+[1 0], [1 1]*Emax, '.b-')
    end
    if ismember(fit_type,[0 1])
        plot(log10(GI50)*[1 1], [0 .5], '.b-')
    end
    
    if p>=.05
        plot(log10(xc), fit_res_flat(xc),'-g');
    end
    
    title(sprintf('r^2 = %.3f', r2))
end


    function [fit_result, gof2] = sigmoidal_fit(doses, response)
        fitopt = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',ranges(1,:),...
            'Upper',ranges(2,:),...
            'Startpoint',priors);
        f = fittype('b + (a-b) / ( 1 + (x/c).^d)','options',fitopt);
        [fit_result,gof2] = fit(doses', response',f);
    end

    function [fit_result, gof2] = flat_fit(doses, response)
        fitopt = fitoptions('Method','NonlinearLeastSquares',...
            'Lower',max(.7,ranges(1,2)),...
            'Upper',ranges(2,1),...
            'Startpoint',priors(1));
        f = fittype('b+0*x','options',fitopt);
        [fit_result,gof2] = fit(doses', response',f);
    end

end