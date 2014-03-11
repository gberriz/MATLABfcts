function [IC50, Hill, Emax, r2, fit_final, p, log] = ...
    ICcurve_fit(Conc, RelativeGrowth, opt)

% parameters : E0 Emax EC50 HS (all in uM, in log10 domain)
priors = [1 .1 1 2];

ranges = [
    .8 1.2  %E0
    0 1    %Emax
    10.^[-4.5 1.5]  %E50
    .1 4    % HS
    ]';

plotting = 0;
fitting = 'average';
fit_type = 'any';

if exist('opt','var')
    fields = {'plotting', 'priors', 'ranges', 'fitting', 'fit_type'};
    for field = fields
        if isfield(opt,field{:})
            eval([field{:} ' = opt.' field{:} ';'])
        end
    end
end

Conc = ToRow(Conc);
assert(all([1 size(RelativeGrowth,2)]==size(Conc)))

% 'IC50'
fit_type = 2;
ranges(1,2) = 0; % lowest Emax


g = mean(RelativeGrowth,1);

[fit_res, gof] = sigmoidal_fit(Conc,g);
[fit_res_flat, gof_flat] = flat_fit(Conc,g);

% F-test for the models
RSS2 = gof.sse;
RSS1 = gof_flat.sse;

df1 = (4 -1);
df2 = (length(g) -4);
F = ( (RSS1-RSS2)/df1 )/( RSS2/df2 );
p = 1-fcdf(F, df1, df2);

xc = 10.^(log10(min(Conc))-1:.05:log10(max(Conc))+1);
r2 = gof.rsquare;

if p>=.05
    GI50 = +Inf;
    IC50 = +Inf;
    Hill = 0;
    Emax = +Inf;
    log = ['**** USING LINEAR FIT **** r2= ' num2str(gof_flat.rsquare,'%.2f')];
    fit_final = fit_res_flat;
    
else
    log = ['r2 = ' num2str(gof.rsquare,'%.2f')];
    
    fit_growth = fit_res(xc);
    
    if any(fit_growth<.5) && any(fit_growth>.5)
        GI50 = fit_res.c*((((fit_res.a-fit_res.b)/(.5-fit_res.b))-1)^(1/fit_res.d));
    elseif all(fit_growth>.5)
        GI50 = Inf;
        log = [log '\tGI50>10*max(Conc) --> +Inf'];
    elseif all(fit_growth<.5)
        GI50 = -Inf;
        log = [log '\tGI50<min(Conc)/10 --> -Inf'];
    else
        warning('undefined GI50')
        log = [log '\tundefined GI50 --> NaN'];
    end
    
    Emax = fit_res.b;
    IC50 = fit_res.c;
    Hill = fit_res.d;
    fit_final = fit_res;
    
end


if plotting
    
    errorbar(log10(Conc), mean(RelativeGrowth,1), std(RelativeGrowth,[],1) ,'.k-');
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
%     if fit_type==1
%         ylim([max(min(min(ylim),-2),-3) max(ylim)])
%     end
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
            'Lower',ranges(1,2),...
            'Upper',ranges(2,1),...
            'Startpoint',priors(1));
        f = fittype('b+0*x','options',fitopt);
        [fit_result,gof2] = fit(doses', response',f);
    end

end