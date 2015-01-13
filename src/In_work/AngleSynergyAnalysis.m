
function [allGIs, t_rfits] = AngleSynergyAnalysis(t_data, extraGIs)
% AngleSynergyAnalysis(t_data)

global Plotting_parameters
%%%%%%%%%%%%%%
%%%% assume outputs from the D300 DrugResponse Suite (like t_mean)
%%%%%%%%%%%%%%
assert(~isempty(t_data),'Empty data')

if isvariable(t_data,'CellLine')
    assert(length(unique(t_data.CellLine))==1, 'There are more than one cell line: %s', ...
        strjoin(cellstr(unique(t_data.CellLine)), ', '))
end

Drugs = unique(t_data.DrugName(t_data.DrugName2~='-'));
assert(~isempty(Drugs), 'No ''Drug1'' found')
assert(length(Drugs)==1, 'There is more than one ''Drug1'': %s', ...
    strjoin(cellstr(Drugs),', '))

temp = setdiff(unique(t_data.DrugName(t_data.DrugName2=='-')), Drugs);
assert(length(temp)==1, 'There is more than one ''Drug2'': %s', ...
    strjoin(cellstr(temp),', '))
Drugs(2) = temp;

%%%%%%%%%%%%%%
%%%% currently assume growth data with 'RelGrowth'
%%%% do NOT handle other metics
%%%%%%%%%%%%%%

%% evaluate the GI50 for each drug independently
%%%%%%  values for GIxx should be specified as an input %%%
if ~exist('extraGIs','var')
    extraGIs = [20 35 65 80 95];
else
    extraGIs = sort(ToRow(extraGIs),'ascend');
end
GIvalues = [50 extraGIs];
GIs = Inf(length(GIvalues),2); % already set to max concentration
Hills = NaN(1,2);
Einfs = Hills;
EC50s = Hills;
Emaxs = Hills;
r2s = Hills;
NormFactor = Hills;

for i=1:2
    t_temp = sortrows(t_data(t_data.DrugName==Drugs(i) & t_data.DrugName2=='-',:), ...
        'Conc');
    Emaxs(i) = 1-min(t_temp.RelGrowth);
    
    opt = struct('pcutoff',.1,'extrapolrange',10^1.5,'Robust',true);
    [GIs(1,i), Hills(i), Einfs(i), ~, r2s(i), EC50s(i), fitfct] = ...
        ICcurve_fit(t_temp.Conc, t_temp.RelGrowth, 'GI50', opt);
    %%%% replace by the actual subfunction and stack it in the file
    %%%%%%%%%%
    
    fprintf('\tfit for %-12s: EC50=%.3g, GI50=%.3g, r2=%.2f\n', char(Drugs(i)), EC50s(i), GIs(1,i), r2s(i));
    warnassert(r2s(i)>.7, 'Bad fit !')
    %%%%% should that stop the algorithm ??
    
    % evaluation of the other GIs
    %%%%% repeated code; could be paste in a subfunction
    fitopts = optimoptions('fsolve', 'Display','none');
    for j=1:length(extraGIs)
        if Einfs(i)>extraGIs(j)/100
            GIs(1+j,i) = fsolve(@(x) fitfct(x)-(1-extraGIs(j)/100), ...
                median(t_temp.Conc), fitopts);
        else
            GIs(1+j,i) = max(t_temp.Conc)*1.5;
        end
    end
    %%%%%% the other GI should directly be evaluated in the ICcurve_fit
    %%%%%% function and capped as the GI50 if the values are extrapolated
    
    
    %%%%%%%%%%%%%%
    %%%%  how to deal with normalization when there is no GI50 for normalization ????
    %%%%  use EC50??? instead (with a warning?)
    %%%%%%%%%%%%%%
    
    if isfinite(GIs(1,1))
        NormConc{i} = 'GI50';
        NormFactor(i) = GIs(1,i);
    else
        NormConc{i} = 'EC50';
        NormFactor(i) = EC50s(i);
        warnprintf('Normalization by EC50 instead of GI50')
    end
    
end


%% caculate the ratios on which to fit the sigmoidal cruves

t_combodata = t_data(t_data.DrugName==Drugs(1) & t_data.DrugName2==Drugs(2),...
    {'Conc' 'Conc2' 'RelGrowth'});

BAratio = 1e-2*round(log10(t_combodata.Conc2./t_combodata.Conc)*1e2);
t_combodata = [t_combodata table(BAratio) ...
    table(t_combodata.Conc+t_combodata.Conc2, 'variablename', {'ConcTot'}) ];
BAratio = unique(t_combodata.BAratio,'sorted');
ratiocounts = hist(t_combodata.BAratio, BAratio);

BAratio = BAratio(ratiocounts>=4);
fprintf('\tFitting for %i ratios:\n', length(BAratio));

t_rfits = [table(BAratio) array2table(NaN(length(BAratio),4), ...
    'variablenames', {'Hill' 'Einf' 'r2' 'Emax'})];
colors = parula(length(BAratio));

ratioGIs = Inf(length(BAratio),length(GIvalues));

get_newfigure(999,[60 400 900 400]);
get_subaxes(1,3,1,[],1)
for i = 1:height(t_rfits)
    t_sub = sortrows(t_combodata(t_combodata.BAratio == t_rfits.BAratio(i),:), ...
        'ConcTot');
    t_rfits.Emax(i) = 1-min(t_sub.RelGrowth);
    
    %%%%%%%
    opt = struct('pcutoff',1,'extrapolrange',10^1.5);
    [ratioGIs(i,1), t_rfits.Hill(i), t_rfits.Einf(i), ~, t_rfits.r2(i), ...
        ~, fitfct, ~, log] = ...
        ICcurve_fit(t_sub.ConcTot, t_sub.RelGrowth, 'GI50', opt);
    %%%% replace by the actual subfunction and stack it in the file
    %%%%%%%%%%
    
    %%%%% repeated code; could be paste in a subfunction
    for j=1:length(extraGIs)
        fitopts = optimoptions('fsolve', 'Display','none');
        temp = fsolve(@(x) fitfct(x)-(1-extraGIs(j)/100), ...
            median(t_sub.ConcTot), fitopts);
        % capping the values to avoid too much extrapolation
        if temp<min(t_sub.ConcTot)*10^-1.5 || imag(temp)~=0
            ratioGIs(i,j+1) = -Inf;
        elseif temp>10^1.5*max(t_sub.ConcTot)
            ratioGIs(i,j+1) = Inf;
        else
            ratioGIs(i,j+1) = temp;
        end
    end
    %%%%%% the other GI should directly be evaluated in the ICcurve_fit
    %%%%%% function and capped as the GI50 if the values are extrapolated
    
    fprintf('\t\tfit for BAratio %-5g: GI50=%.3g, r2=%.2f\n', t_rfits.BAratio(i), ...
        ratioGIs(i,1), t_rfits.r2(i));
    
    plot(t_sub.ConcTot, t_sub.RelGrowth, 'o', 'color', colors(i,:))
    
    if isinf(ratioGIs(i,1)) || imag(ratioGIs(i,1))~=0
        if imag(ratioGIs(i,1))~=0
            ratioGIs(i,1) = -Inf;
        end
        fprintf(['\t\t --> ' log '\n'])
    end
    x = 10.^(log10(infmin([ratioGIs(i,1);t_sub.ConcTot])):.01:...
        log10(infmax([ratioGIs(i,1);t_sub.ConcTot])));
    plot(x, fitfct(x), '-', 'color', colors(i,:))
    plot(x(end)*[1 1.2], [1 1]-t_rfits.Emax(i), ':', 'color', colors(i,:))
    plot(x(end)*[1 1.2], [1 1]-t_rfits.Einf(i), '--', 'color', colors(i,:))
    
    plot(ratioGIs(i,:), 1-GIvalues/100, '.', 'color', colors(i,:), ...
        'markersize', 15)
end
set(gca,'xscale','log')
t_rfits.RelBAratio = t_rfits.BAratio +log10(NormFactor(1)) -log10(NormFactor(2));

%% calculate the GIs along the directions of sampling (no Hill or Emax here)
Conc = cell(2,1);
AlignConc = cell(2,1);
for i=1:2
    Conc{i} = unique(t_data.Conc(t_data.DrugName==Drugs(i)),'sorted');
end

alignedGIs = cell(2,1);

for iC=1:2
    
    if iC==1 % not the most elegant solution
        n = hist(t_combodata.Conc,Conc{iC});
    else
        n = hist(t_combodata.Conc2,Conc{iC});
    end
    idx = find(n>=4);
    AlignConc{iC} = Conc{iC}(idx);
    alignedGIs{iC} = Inf(length(idx), length(GIvalues));
    
    colors = parula(length(idx));
    get_subaxes(1,3,1+iC,[],1)
    for i=1:length(idx)
        
        if iC==1 % swapping to keep code the same
            t_sub = t_combodata(t_combodata.Conc==Conc{iC}(idx(i)),:);
        else
            t_sub = t_combodata(t_combodata.Conc2==Conc{iC}(idx(i)),:);
            t_sub.Properties.VariableNames({'Conc' 'Conc2'}) = {'Conc2' 'Conc'};
        end
        
        RelGrowth0 = mean(t_data.RelGrowth(t_data.DrugName==Drugs(iC) & t_data.DrugName2=='-' & ...
            t_data.Conc==Conc{iC}(idx(i))));
        
        ranges = [
            RelGrowth0-.05 max([RelGrowth0;t_sub.RelGrowth])+.05  %E0
            0 1    %Emax
            max(min(t_sub.Conc2)*1e-4,1e-7) min(max(t_sub.Conc2)*1e2, 1e3)  %E50
            .1 5    % HS
            ]';
        opt = struct('pcutoff',1,'extrapolrange',10^1.5,'ranges',ranges, ...
            'Robust', true);
        [alignedGIs{iC}(i,1), ~, ~, ~, r2, ~, fitfct, ~, log] = ...
            ICcurve_fit(t_sub.Conc2, t_sub.RelGrowth, 'GI50', opt);
        
        if r2<.75
            fprintf(['\t\t Bad fit --> ' log '\n'])
            alignedGIs{iC}(i,:) = NaN;
            continue
        end
        
        %%%%% repeated code; could be paste in a subfunction
        for j=1:length(extraGIs)
            fitopts = optimoptions('fsolve', 'Display','none');
            [temp, fval, flag] = fsolve(@(x) fitfct(x)-(1-extraGIs(j)/100), ...
                median(t_sub.Conc2), fitopts);
            if flag>=0 % convergence
                % capping the values to avoid too much extrapolation
                if temp<min(t_sub.Conc2)*10^-1.5
                    alignedGIs{iC}(i,j+1) = -Inf;
                elseif temp>10^1.5*max(t_sub.Conc2)
                    alignedGIs{iC}(i,j+1) = Inf;
                else
                    alignedGIs{iC}(i,j+1) = temp;
                end
            elseif min(t_sub.RelGrowth)>(1-extraGIs(j)/100)
                    alignedGIs{iC}(i,j+1) = Inf;
            elseif max(t_sub.RelGrowth)<(1-extraGIs(j)/100)
                    alignedGIs{iC}(i,j+1) = -Inf;
            else
                error('unknown case for fsolve')
            end
                
                
        end
        
        
        fprintf('\t\tfit for Conc%i %-5g: GI50=%.3g, r2=%.2f\n', iC, Conc{iC}(idx(i)), ...
            alignedGIs{iC}(i,1), r2);
        
        plot(t_sub.Conc2, t_sub.RelGrowth, 'o', 'color', colors(i,:))
        
        if isinf(alignedGIs{iC}(i,1)) || imag(alignedGIs{iC}(i,1))~=0
            if imag(alignedGIs{iC}(i,1))~=0
                alignedGIs{iC}(i,1) = -Inf;
            end
            fprintf(['\t\t --> ' log '\n'])
        end
        x = 10.^(log10(infmin([alignedGIs{iC}(i,1);t_sub.Conc2])):.01:...
            log10(infmax([alignedGIs{iC}(i,1);t_sub.Conc2])));
        
        plot(x, fitfct(x), '-', 'color', colors(i,:))        
        plot(x(1)*[1 1.2], [1 1]*mean(ranges(:,1)), '--', 'color', colors(i,:))
        
        plot(alignedGIs{iC}(i,:), 1-GIvalues/100, '.', 'color', colors(i,:), ...
            'markersize', 15)
        
    end
    set(gca,'xscale','log')
    
end



%% getting the actual Concentration for the ratioGIxx for plotting later
% and stacking with the alignedGIs

allGIs = struct('GIval', num2cell(GIvalues), 'BAratio', [], 'RealConc', [], ...
    'NullConc', [], 'CI', [], 'RelBAratio', []);

w = [ones(height(t_rfits),1) 10.^t_rfits.BAratio];
w = w ./ (sum(w,2)*[1 1]);


% smoothing parameter
Nextra = 4;
Nsmooth = 5;
    
get_newfigure(998,[750 50 900 700])
for i=1:length(GIvalues)
    %%
    BAratio = real([ t_rfits.BAratio;
        log10(alignedGIs{1}(:,i)./AlignConc{1});
        log10(AlignConc{2}./alignedGIs{2}(:,i))]);
    RealConc = [repmat(ratioGIs(:,i),1,2).*w; 
        [AlignConc{1} alignedGIs{1}(:,i)] 
        [alignedGIs{2}(:,i) AlignConc{2}]];
    
    idx = find(isfinite(BAratio) & ...
        BAratio<log10((max(t_combodata.Conc2)/min(t_combodata.Conc))) & ...
        BAratio>log10((min(t_combodata.Conc2)/max(t_combodata.Conc))) & ...
        all(isfinite(RealConc),2));
        
    idx = idx(sortidx(BAratio(idx)));
    allGIs(i).BAratio = BAratio(idx)';
    
    allGIs(i).RealConc = RealConc(idx,:);
    w2 = [ones(length(allGIs(i).BAratio),1) 10.^allGIs(i).BAratio'];
    w2 = w2 ./ (sum(w2,2)*[1 1]);    
    
    % allGIs(i).CI = sum(allGIs(i).RealConc,2) ./ (w2*GIs(i,:)');
    allGIs(i).CI = sum((allGIs(i).RealConc ./ repmat(GIs(i,:), length(allGIs(i).RealConc),1)),2) ./ ...
        sum((w2.*repmat(GIs(i,:),length(w2),1))./repmat(GIs(i,:), length(allGIs(i).RealConc),1),2);
    
    NullBAratio = ((infmin(BAratio)-1.5):.1:(infmax(BAratio)+1.5));
    w2 = [ones(length(NullBAratio),1) 10.^NullBAratio'];
    w2 = w2 ./ (sum(w2,2)*[1 1]);    
    allGIs(i).NullConc = w2.*repmat(GIs(i,:),length(w2),1);
    
    
    
    % smoothing
    Nextra = 4;
    BAratio = [.1*floor(min(allGIs(i).BAratio)*10)-((Nextra-1):-1:0)*.1-.5 ...
        allGIs(i).BAratio ...
        .1*ceil(max(allGIs(i).BAratio)*10)+(0:(Nextra-1))*.1+.5 ];
    allGIs(i).intrpBAratio = ...
        ((infmin(allGIs(i).BAratio)-.25):.2:(infmax(allGIs(i).BAratio)+.25))';
    
    
    c1 = [repmat(GIs(i,1),Nextra,1);  ...
        allGIs(i).RealConc(:,1); ...
        -Inf(Nextra ,1)];
    c1(c1==-Inf) = infmin([t_combodata.Conc;allGIs(i).RealConc(:,1)])*10^-1.5;
%     c1(c1==Inf) = infmax([t_combodata.Conc;allGIs(i).RealConc(:,1)])*10^1.5;
    c2 = [-Inf(Nextra ,1); ...
        allGIs(i).RealConc(:,2); ...
        repmat(GIs(i,2),Nextra ,1)];
    c2(c2==-Inf) = infmin([t_combodata.Conc2;allGIs(i).RealConc(:,2)])*10^-1.5;
%     c2(c2==Inf) = infmax([t_combodata.Conc2;allGIs(i).RealConc(:,2)])*10^1.5;
    
    allGIs(i).intrpConc = 10.^[smooth(interp1(BAratio, log10(c1), ...
        allGIs(i).intrpBAratio, 'PCHIP'),Nsmooth) ...
        smooth(interp1(BAratio, log10(c2), ...
        allGIs(i).intrpBAratio, 'PCHIP'),Nsmooth)];
    
    CI = [ones(Nextra ,1); allGIs(i).CI; ones(Nextra ,1)];
    allGIs(i).intrpCI = 2.^smooth(interp1(BAratio, log2(CI), ...
        allGIs(i).intrpBAratio, 'PCHIP'),Nsmooth);
    
    % plot QC
    get_subaxes(length(GIvalues),3,i,1,1)
    plot(allGIs(i).BAratio, log10(allGIs(i).RealConc(:,1)), 'x')
    plot(BAratio, log10(c1), '.', 'markersize',12)
    plot(allGIs(i).intrpBAratio, log10(allGIs(i).intrpConc(:,1)), '-')
    plot(NullBAratio, log10(allGIs(i).NullConc(:,1)), '--k')
    plot(-log10(NormFactor(1))+log10(NormFactor(2))*[1 1], [-2 1], '-k')
    plot([0 0], [-2 1], '-', 'color', [.5 .5 .5])
    ylim([-3.5 1.5])
    
    get_subaxes(length(GIvalues),3,i,2,1)
    plot(allGIs(i).BAratio, log10(allGIs(i).RealConc(:,2)), 'x')
    plot(BAratio, log10(c2), '.', 'markersize',12)
    plot(allGIs(i).intrpBAratio, log10(allGIs(i).intrpConc(:,2)), '-')
    plot(NullBAratio, log10(allGIs(i).NullConc(:,2)), '--k')
    plot(-log10(NormFactor(1))+log10(NormFactor(2))*[1 1], [-2 1], '-k')
    plot([0 0], [-2 1], '-', 'color', [.5 .5 .5])
    ylim([-3.5 1.5])
    
    
    get_subaxes(length(GIvalues),3,i,3,1)
    plot(allGIs(i).BAratio, log2(allGIs(i).CI), 'x')
    plot(BAratio, log2(CI), '.', 'markersize',12)
    plot(allGIs(i).intrpBAratio, log2(allGIs(i).intrpCI), '-')
    plot(-log10(NormFactor(1))+log10(NormFactor(2))*[1 1], [-2 1], '-k')
    plot([0 0], [-2 1], '-', 'color', [.5 .5 .5])
    ylim([-3.5 1.5])
   
    
    allGIs(i).RelBAratio = allGIs(i).BAratio +log10(NormFactor(1)) -log10(NormFactor(2));
    allGIs(i).intrpRelBAratio = allGIs(i).intrpBAratio +log10(NormFactor(1)) -log10(NormFactor(2));
end

%% smoothing of Hill coeff and Emax/inf
intrpxlims = [.1*floor(min([allGIs.RelBAratio])*10-1) ...
    .1*ceil(max([allGIs.RelBAratio])*10+1)];

RelBAratio = [intrpxlims(1)-(1:(Nextra-1))*.2 ...
    t_rfits.RelBAratio' intrpxlims(2)+(1:(Nextra-1))*.2];
intrpRelrat = (intrpxlims(1)-.25):.2:(intrpxlims(2)+.25);
intrpHill = smooth(interp1(RelBAratio, ...
    [repmat(Hills(1),1,Nextra) t_rfits.Hill' repmat(Hills(2),1,Nextra)], ...
        intrpRelrat, 'PCHIP'),Nsmooth);
intrpEinf = smooth(interp1(RelBAratio, ...
    [repmat(Einfs(1),1,Nextra) t_rfits.Einf' repmat(Einfs(2),1,Nextra)], ...
        intrpRelrat, 'PCHIP'),Nsmooth);
intrpEmax = smooth(interp1(RelBAratio, ...
    [repmat(Emaxs(1),1,Nextra) t_rfits.Emax' repmat(Emaxs(2),1,Nextra)], ...
        intrpRelrat, 'PCHIP'),Nsmooth);



%%

if isvariable(t_data,'CellLine')
get_newfigure(100,[50 50 700 550], [char(unique(t_data.CellLine)) ...
    '_Combo_' char(Drugs(1)) '+' char(Drugs(2)) '.pdf'])
else
get_newfigure(100,[50 50 700 550], ...
    [ 'Combo_' char(Drugs(1)) '+' char(Drugs(2)) '.pdf'])
end
    

% spacing for checkboard
xwidth = .02;
ywidth = .03;
space = .005;

% position of axes
pos = [
    .07 .59 .5 .34; % CI plot
    .67 .59 .3 .34; % Hill plot
    .67 .12 .3 .34]; % Emax plot
pos(4,:) = [pos(1,1) pos(3,2)-ywidth-space .4 pos(3,4)+ywidth]; % checkboard

% position for checkboard
subpos = [
    pos(4,1) pos(4,2)+ywidth+space xwidth pos(4,4)-ywidth;  % Drug 1
    pos(4,1)+xwidth+space pos(4,2) pos(4,3)-xwidth ywidth;  % Drug 2
    pos(4,1)+xwidth+space pos(4,2)+ywidth+space pos(4,3)-xwidth pos(4,4)-ywidth; %combo
    pos(4,1)+pos(4,3)+5*space pos(4,2)+.6*pos(4,4) xwidth/2 pos(4,4)*.4; % colorbar
    pos(4,1)+pos(4,3)+2*space pos(4,2) 0 0]; % legend

% define the ticks for the ratio axis
xticks = [intrpxlims(1)-.1 (intrpxlims(1)+.2):(.1*round(diff(intrpxlims)*1.2)):(intrpxlims(2)-.2) ...
    intrpxlims(2)+.1];
xlims = [intrpxlims(1)-.1 intrpxlims(2)+.1];
Concxticks = xticks -log10(NormFactor(1)) +log10(NormFactor(2));
% ticks in ratio of concentration
Concxticklabels = cell(length(xticks),1);
for i=2:length(xticks)-1
    if Concxticks(i)<0
        Concxticklabels{i} = sprintf('%.1f:1', 10^(-Concxticks(i)));
    else
        Concxticklabels{i} = sprintf('1:%.1f', 10^(Concxticks(i)));
    end
end
Concxticklabels{1} = char(Drugs(1));
Concxticklabels{end} = char(Drugs(2));

% ticks in ratio relative to GI50/EC50
Relxticklabels = cell(length(xticks),1);
for i=2:length(xticks)-1
    if xticks(i)<0
        Relxticklabels{i} = sprintf('%.1f:1', 10^(-xticks(i)));
    else
        Relxticklabels{i} = sprintf('1:%.1f', 10^(xticks(i)));
    end
end
Relxticklabels{1} = char(Drugs(1));
Relxticklabels{end} = char(Drugs(2));
xvals = [xticks(1) t_rfits.RelBAratio' xticks(end)];

GIcolors = min([0 0 0; .2+(sum(extraGIs<50):-1:1)'*[0 .2 .8]/sum(extraGIs<50);
    .2+(1:sum(extraGIs>50))'*[.4 .8 0]/sum(extraGIs>50)],1);



%%%%%%%%%%%%%%%
% first plot for the CI
% conc ratio axis
get_newaxes(pos(1,:),1)
set(gca,'xtick',xticks,'xticklabel',Concxticklabels,'xticklabelrotation',90, ...
    'ytick',[],'xaxislocation','top',Plotting_parameters.axes{:})
xlim(xlims)
xlabel('Drug ratios (Conc)',Plotting_parameters.axislabel{:})


get_newaxes(pos(1,:),1)
plot(repmat([xvals([1 end]) NaN],1,3), reshape([1;1;NaN]*[-1 0 1],[],1)', ...
    '-', 'color', [.7 .7 .7])
plot([0 0], [-32 32], '-k')
plot(log10(NormFactor(1))-log10(NormFactor(2))*[1 1], [-32 32], '-', 'color', [.5 .5 .5])

h = NaN(length(GIvalues),1);
for j=[2:length(GIvalues) 1] % plot GI50 last
    h(j) = plot(allGIs(j).RelBAratio, log2(allGIs(j).CI), ...
        '.','color',GIcolors(j,:), 'markersize', 12);
    h(j) = plot(allGIs(j).intrpRelBAratio, log2(allGIs(j).intrpCI), ...
        '-','color',GIcolors(j,:),'linewidth',1);
end

set(gca,'xtick',xticks,'xticklabel',Relxticklabels,'xticklabelrotation',90, ...
    'ytick',-5:5,'yticklabel',[strcat('1/',num2cellstr(2.^(5:-1:1))) ...
    strcat(num2cellstr(2.^(0:5)),'/1')],Plotting_parameters.axes{:}, 'box','on')
ylim(log2([min([.45 cellfun(@min,{allGIs.CI})]) max([2.2 cellfun(@max,{allGIs.CI})])]))
xlim(xlims)
ylabel('Combination index',Plotting_parameters.axislabel{:})
xlabel(['Drug ratios (' strjoin(NormConc, '/') ' normalized)'],...
    Plotting_parameters.axislabel{:})
% legend(h,[strcat('GI', num2cellstr(GIvalues)) 'Failed points'],'location','best', ...
%     Plotting_parameters.axislabel{:})


%%%%%%%%%%%%%%%%%%%
% plot for the Hill coeff
% get_newaxes(pos(2,:),1)
% set(gca,'xtick',xticks,'xticklabel',Concxticklabels,'xticklabelrotation',90, ...
%     'ytick',[],'xaxislocation','top',Plotting_parameters.axes{:})
% xlim(xlims)
get_newaxes(pos(2,:),1)
plot([0 0], [-32 32], '-k')
plot(log10(NormFactor(1))-log10(NormFactor(2))*[1 1], [-32 32], '-', 'color', [.5 .5 .5])
plot(t_rfits.RelBAratio, t_rfits.Hill, '.k')
plot(intrpRelrat, intrpHill , '-k')
set(gca,'xtick',xticks,'xticklabel',Relxticklabels,'xticklabelrotation',90,...
    Plotting_parameters.axes{:}, 'box','on')
xlim(xlims)
ylabel('Hill coefficient',Plotting_parameters.axislabel{:})
xlabel(['Drug ratios (' strjoin(NormConc, '/') ' normalized)'],...
    Plotting_parameters.axislabel{:})
ylim([min([.4 t_rfits.Hill']) max([4 t_rfits.Hill'])])


%%%%%%%%%%%%%%%%%%%%%%%
% plot for the Emax/inf
% get_newaxes(pos(3,:),1)
% set(gca,'xtick',xticks,'xticklabel',Concxticklabels,'xticklabelrotation',90, ...
%     'ytick',[],'xaxislocation','top',Plotting_parameters.axes{:})
% xlim(xlims)
get_newaxes(pos(3,:),1)
plot([0 0], [-32 32], '-k')
plot(log10(NormFactor(1))-log10(NormFactor(2))*[1 1], [-32 32], '-', 'color', [.5 .5 .5])
plot(xvals([1 end]), [0 0], '-', 'color', [.7 .7 .7])
h = plot(xvals, 1-[Einfs(1) t_rfits.Einf' Einfs(2)], 'xk');
h(2) = plot(xvals, 1-[Emaxs(1) t_rfits.Emax' Emaxs(2)], '.k');
h = plot(intrpRelrat, 1-intrpEinf, '--k');
h(2) = plot(intrpRelrat, 1-intrpEmax, '-k');
set(gca,'xtick',xticks,'xticklabel',Relxticklabels,'xticklabelrotation',90,...
    Plotting_parameters.axes{:}, 'box','on')
xlim(xlims)
ylim([min([-.3 1-t_rfits.Emax' 1-Emaxs]) max([.3 1-t_rfits.Emax' 1-Emaxs])])
ylabel('Lowest growth (1-E_{inf}/E_{max})',Plotting_parameters.axislabel{:})
xlabel(['Drug ratios (' strjoin(NormConc, '/') ' normalized)'],...
    Plotting_parameters.axislabel{:})
legend(h,{'1-E_{inf}' '1-E_{max}'},'location','best',...
    Plotting_parameters.axislabel{:})



% plotting the checkboard for combo

%%%%%%% put the following in its own function

%%%%%% assuming carthesian product on conc; need to be generalized
% needs to be split in the columns for single dose and matrix for combos
% (for cases where there is combo are not carthesian product on conc)
%%%
logC1 = log10(Conc{1});
logC2 = log10(Conc{2});
GrowthMx = NaN(length(logC1)+1, length(logC2)+1);
GrowthMx(1,1) = 1;

for iC1 = 1:length(logC1)
    GrowthMx(iC1+1,1) = t_data.RelGrowth(t_data.DrugName==Drugs(1) & ...
        t_data.Conc==Conc{1}(iC1) & t_data.DrugName2=='-');
    
    for iC2 = 1:length(logC2)
        if iC1==1
            GrowthMx(1,iC2+1) = t_data.RelGrowth(t_data.DrugName==Drugs(2) & ...
                t_data.Conc==Conc{2}(iC2) & t_data.DrugName2=='-');
        end
        
        val = t_data.RelGrowth(t_data.DrugName==Drugs(1) & ...
            t_data.Conc==Conc{1}(iC1) & t_data.DrugName2==Drugs(2) & ...
            t_data.Conc2==Conc{2}(iC2) );
        if ~isempty(val)
            GrowthMx(iC1+1,iC2+1) = val;
        end
        
        
    end
end
%%%%%%
% the combo should be an interpolation (like pcolor) in the log 10
% (for cases where the doses are not equidistant in the log domain)
%%%%%%


xlims = [1.5*logC2(1)-.5*logC2(2) 1.5*logC2(end)-.5*logC2(end-1)];
ylims = [1.5*logC1(1)-.5*logC1(2) 1.5*logC1(end)-.5*logC1(end-1)];


% small corner (control)
get_newaxes([pos(4,1) pos(4,2) xwidth ywidth])
imagesc(1,1,1,[-1 1])
set(gca,'xtick',1,'xticklabel',0, 'ytick',1,'yticklabel',0 , ...
    Plotting_parameters.axes{:},'ydir','normal')

% plot for Drug 1
get_newaxes(subpos(1,:),1)
imagesc(ones(1,length(logC1)), logC1, GrowthMx(2:end,1),[-1 1]);
ylabel([char(Drugs(1))  '[uM]'], Plotting_parameters.axislabel{:})
set(gca,'xtick',[], 'ytick', logC1(2:2:end), ...
    'yticklabel', num2cellstr(Conc{1}(2:2:end),'%.2g'), ...
    Plotting_parameters.axes{:},'ydir','normal')
% plot the GIxx for Drug 1
for i=1:length(GIvalues)
    plot([.5 1.5], log10(GIs(i,1))*[1 1],'-','color', GIcolors(i,:), ...
        'linewidth',2);
end
ylim(ylims)

% plot for Drug 2
get_newaxes(subpos(2,:),1)
imagesc(logC2, ones(1,length(logC2)), GrowthMx(1,2:end),[-1 1]);
xlabel([char(Drugs(2)) '[uM]'], Plotting_parameters.axislabel{:})
set(gca,'ytick',[], 'xtick', logC2(2:2:end), ...
    'xticklabel', num2cellstr(Conc{2}(2:2:end),'%.2g'),...
    Plotting_parameters.axes{:},'ydir','normal')
% plot the GIxx for Drug 2
for i=1:length(GIvalues)
    plot(log10(GIs(i,2))*[1 1], [.5 1.5], '-','color', GIcolors(i,:), ...
        'linewidth',2);
end
xlim(xlims)

% plot for combination
get_newaxes(subpos(3,:),1)
imagesc(logC2, logC1, GrowthMx(2:end,2:end),[-1 1]);
set(gca,'xtick',[], 'ytick',[],'ydir','normal')
% set(gca, 'ytick', logC1(2:2:end), ...
%     'yticklabel', num2cellstr(Conc{1}(2:2:end),'%.2g'), ...
%      'xtick', logC2(2:2:end),'xticklabel', num2cellstr(Conc{2}(2:2:end),'%.2g'),...
%     Plotting_parameters.axes{:},'ydir','normal')

% plot the ratio directions used for fit
for i=1:height(t_rfits)
    plot([-5 5], [-5 5]-t_rfits.BAratio(i), '-', 'color', [.8 .8 .8]);
end

% plot the isoGI cruves
h = NaN(length(GIvalues)+1,1);
for j=[ 2:length(GIvalues) 1]    
    % theoric curve (additive)
    h(j) = plot(log10(allGIs(j).NullConc(:,2)), log10(allGIs(j).NullConc(:,1)), ...
        '--','color', GIcolors(j,:), 'linewidth',1);
    % real curve
    h(j) = plot(log10(allGIs(j).RealConc(:,2)), log10(allGIs(j).RealConc(:,1)), ...
        '.','color', GIcolors(j,:), 'markersize', 12);
    h(j) = plot(log10(allGIs(j).intrpConc(:,2)), log10(allGIs(j).intrpConc(:,1)), ...
        '-','color', GIcolors(j,:), 'linewidth',1);
    
end
h(end) = plot([NaN NaN], [NaN NaN], '--', 'color', [.3 .3 .3]);


% plot the ratio directions (1:1) and (GI50:GI50)
plot([min([logC1;logC2])-5 max([logC1;logC2])+5], ...
    [min([logC1;logC2])-5 max([logC1;logC2])+5], '-', 'color', [.5 .5 .5]);

plot([min([logC1;logC2])-5 max([logC1;logC2])+5]-log10(NormFactor(1)), ...
    [min([logC1;logC2])-5 max([logC1;logC2])+5]-log10(NormFactor(2)), '-k');


xlim(xlims)
ylim(ylims)
colormap(Plotting_parameters.cmapGrWBr)

hc = colorbar('location','east');
set(hc,'position',subpos(4,:))
ylabel(hc,'RelGrowth', Plotting_parameters.axislabel{:})
ylim(hc,[-.3 1])

hl = legend(h, [strcat('GI', num2cellstr(GIvalues)) 'CI=1'], ...
    Plotting_parameters.axislabel{:});
posl = get(hl,'position');
set(hl,'position', [subpos(5,1) subpos(5,2) posl(3:4)])