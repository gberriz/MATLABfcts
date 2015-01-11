
function AngleSynergyAnalysis(t_data)
% AngleSynergyAnalysis(t_data)

global Plotting_parameters
%%%%%%%%%%%%%%
%%%% assume outputs from the D300 DrugResponse Suite (like t_mean)
%%%%%%%%%%%%%%
assert(~isempty(t_data),'Empty data')

assert(length(unique(t_data.CellLine))==1, 'There are more than one cell line: %s', ...
    strjoin(cellstr(unique(t_data.CellLine)), ', '))

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
extraGIs = sort(ToRow([20 35 65 80 95]),'ascend');
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
    
    opt = struct('pcutoff',.1,'extrapolrange',30,'Robust',true);
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

BAratio = BAratio(ratiocounts>4);
fprintf('\tFitting for %i ratios:\n', length(BAratio));

t_rfits = [table(BAratio) array2table(NaN(length(BAratio),4), ...
    'variablenames', {'Hill' 'Einf' 'r2' 'Emax'})];
colors = parula(length(BAratio));

ratioGIs = Inf(length(BAratio),length(GIvalues));

get_newfigure(999);
get_subaxes(1,3,1,[],1)
for i = 1:height(t_rfits)
    t_sub = sortrows(t_combodata(t_combodata.BAratio == t_rfits.BAratio(i),:), ...
        'ConcTot');
    t_rfits.Emax(i) = 1-min(t_sub.RelGrowth);
    
    %%%%%%%
    opt = struct('pcutoff',1,'extrapolrange',30);
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
        if temp<1e-2*min(t_sub.ConcTot)
            ratioGIs(i,j+1) = -Inf;
        elseif temp>1e2*max(t_sub.ConcTot)
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
    
    if isinf(ratioGIs(i,1))
        fprintf(['\t\t --> ' log '\n'])
    else
        x = 10.^(log10(min([ratioGIs(i,1);t_sub.ConcTot])):.01:...
            log10(max([ratioGIs(i,1);t_sub.ConcTot])));
        plot(x, fitfct(x), '-', 'color', colors(i,:))
        plot(x(end)*[1 1.2], [1 1]-t_rfits.Emax(i), ':', 'color', colors(i,:))
        plot(x(end)*[1 1.2], [1 1]-t_rfits.Einf(i), '--', 'color', colors(i,:))
    end
end
set(gca,'xscale','log')


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
    idx = find(n>4);
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
        
        ranges = [
            t_data.RelGrowth(t_data.DrugName==Drugs(iC) & t_data.DrugName2=='-' & ...
            t_data.Conc==Conc{iC}(idx(i)))+[-.05 .05]  %E0
            0 1    %Emax
            max(min(t_sub.Conc2)*1e-4,1e-7) min(max(t_sub.Conc2)*1e2, 1e3)  %E50
            .1 5    % HS
            ]';
        opt = struct('pcutoff',1,'extrapolrange',30,'ranges',ranges);
        [alignedGIs{iC}(i,1), ~, ~, ~, r2, ~, fitfct, ~, log] = ...
            ICcurve_fit(t_sub.Conc2, t_sub.RelGrowth, 'GI50', opt);
        
        %%%%% repeated code; could be paste in a subfunction
        for j=1:length(extraGIs)
            fitopts = optimoptions('fsolve', 'Display','none');
            temp = fsolve(@(x) fitfct(x)-(1-extraGIs(j)/100), ...
                median(t_sub.ConcTot), fitopts);
            % capping the values to avoid too much extrapolation
            if temp<1e-2*min(t_sub.ConcTot) || imag(temp)~=0
                alignedGIs{iC}(i,j+1) = -Inf;
            elseif temp>1e2*max(t_sub.ConcTot)
                alignedGIs{iC}(i,j+1) = Inf;
            else
                alignedGIs{iC}(i,j+1) = temp;
            end
        end
        
        
        fprintf('\t\tfit for Conc%i %-5g: GI50=%.3g, r2=%.2f\n', iC, Conc{iC}(idx(i)), ...
            alignedGIs{iC}(i,1), r2);
        
        plot(t_sub.Conc2, t_sub.RelGrowth, 'o', 'color', colors(i,:))
        
        if isinf(alignedGIs{iC}(i,1)) || imag(alignedGIs{iC}(i,1))~=0
            alignedGIs{iC}(i,1) = -Inf;
            fprintf(['\t\t --> ' log '\n'])
        else
            x = 10.^(log10(min([alignedGIs{iC}(i,1);t_sub.Conc2])):.01:...
                log10(max([alignedGIs{iC}(i,1);t_sub.Conc2])));
            plot(x, fitfct(x), '-', 'color', colors(i,:))
            
            plot(x(1)*[1 1.2], [1 1]*mean(ranges(:,1)), '--', 'color', colors(i,:))
        end
    end
    set(gca,'xscale','log')
    
end

% for i=1:2
%     GIConc{i} = min(GIConc{i}, max(t_data.Conc(t_data.DrugName==Drugs(i)))*100);
% end


%% getting the actual Concentration for the ratioGIxx for plotting later
%% and stacking with the alignedGIs
w = [ones(height(t_rfits),1) 10.^t_rfits.BAratio];
w = w ./ (sum(w,2)*[1 1]);

GIConc = cell(2,1);
NullGIConc = GIConc;
for i=1:2
    %GIConc{i} = w(:,i)*ratioGIs(:,i)';
    GIConc{i} = ratioGIs.*repmat(w(:,i),1,length(GIvalues));
    NullGIConc{i} = w(:,i)*GIs(:,i)';
end


allGIs = struct('GIval', num2cell(GIvalues), 'BAratio', [], 'RealConc', [], ...
    'NullConc', [], 'CI', [], 'RelBAratio', []);

for i=1:length(GIvalues)
    BAratio = [ t_rfits.BAratio;
        log10(alignedGIs{1}(:,i)./AlignConc{1});
        log10(AlignConc{2}./alignedGIs{2}(:,i))];
    RealConc = [repmat(ratioGIs(:,i),1,2).*w; 
        [AlignConc{1} alignedGIs{1}(:,i)] 
        [alignedGIs{2}(:,i) AlignConc{2}]];
    w2 = [ones(length(BAratio),1) 10.^BAratio];
    w2 = w2 ./ (sum(w2,2)*[1 1]);
    
    NullConc = w2.*repmat(GIs(i,:),length(w2),1);
    
    CI = [sum(RealConc,2) ./ (w2*GIs(i,:)')];
    
    idx = find(isfinite(BAratio) & imag(BAratio)==0);
    idx = idx(sortidx(BAratio(idx)));
    allGIs(i).BAratio = BAratio(idx);
    allGIs(i).RelBAratio = allGIs(i).BAratio +log10(NormFactor(1)) -log10(NormFactor(2));
    allGIs(i).RealConc = RealConc(idx,:);
    allGIs(i).NullConc = NullConc(idx,:);
    allGIs(i).CI = CI(idx);
end


% 
% GIConc{1} = [GIs(:,1)'; GIConc{1}; ...
%     repmat(min(t_data.Conc(t_data.DrugName==Drugs(2)))/100,1,length(GIvalues))];
% GIConc{2} = [repmat(min(t_data.Conc(t_data.DrugName==Drugs(1)))/100,1,length(GIvalues)); ...
%     GIConc{2}; GIs(:,2)'];
% 
% NullGIConc{1} = [GIs(:,1)'; NullGIConc{1}; ...
%     repmat(min(t_data.Conc(t_data.DrugName==Drugs(2)))/100,1,length(GIvalues))];
% NullGIConc{2} = [repmat(min(t_data.Conc(t_data.DrugName==Drugs(1)))/100,1,length(GIvalues)); ...
%     NullGIConc{2}; GIs(:,2)'];
% 
% 
% 
% %% evaluate the CI
CI = ratioGIs ./ (w*GIs');

%%%%%%%%%%%%%%
%%%%  how to deal with normalization when there is no GI50 for normalization ????
%%%%  use EC50??? instead (with a warning?)
%%%%%%%%%%%%%%

t_rfits.RelBAratio = t_rfits.BAratio +log10(NormFactor(1)) -log10(NormFactor(2));


%%

get_newfigure(100,[50 50 700 550], [char(unique(t_data.CellLine)) ...
    '_Combo_' char(Drugs(1)) '+' char(Drugs(2)) '.pdf'])

% spacing for checkboard
xwidth = .03;
ywidth = .04;
space = .005;

% position of axes
pos = [
    .07 .6 .5 .37; % CI plot
    .67 .6 .3 .37; % Hill plot
    .67 .12 .3 .37]; % Emax plot
pos(4,:) = [pos(1,1) pos(3,2)-ywidth-space .4 pos(3,4)+ywidth]; % checkboard

% position for checkboard
subpos = [
    pos(4,1) pos(4,2)+ywidth+space xwidth pos(4,4)-ywidth;  % Drug 1
    pos(4,1)+xwidth+space pos(4,2) pos(4,3)-xwidth ywidth;  % Drug 2
    pos(4,1)+xwidth+space pos(4,2)+ywidth+space pos(4,3)-xwidth pos(4,4)-ywidth; %combo
    pos(4,1)+pos(4,3)+5*space pos(4,2)+.5*pos(4,4) xwidth/2 pos(4,4)/2; % colorbar
    pos(4,1)+pos(4,3)+2*space pos(4,2) 0 0]; % legend



% define the ticks for the ratio axis
if height(t_rfits)>5
    tickidx = 2:floor(height(t_rfits)/5):(height(t_rfits)-1);
else
    tickidx = height(t_rfits);
end

% ticks in ratio of concentration
xticks = [2*t_rfits.BAratio(1)-t_rfits.BAratio(2)
    t_rfits.BAratio(tickidx)
    2*t_rfits.BAratio(end)-t_rfits.BAratio(end-1)];
xticklabels = cell(length(xticks),1);
for i=2:length(xticks)-1
    if xticks(i)<0
        xticklabels{i} = sprintf('%.1f:1', 10^(-xticks(i)));
    else
        xticklabels{i} = sprintf('1:%.1f', 10^(xticks(i)));
    end
end
xticklabels{1} = char(Drugs(1));
xticklabels{end} = char(Drugs(2));

% ticks in ratio relative to GI50/EC50
Relxticks = [2*t_rfits.RelBAratio(1)-t_rfits.RelBAratio(2)
    t_rfits.RelBAratio(tickidx)
    2*t_rfits.RelBAratio(end)-t_rfits.RelBAratio(end-1)];
Relxticklabels = cell(length(Relxticks),1);
for i=2:length(Relxticks)-1
    if Relxticks(i)<0
        Relxticklabels{i} = sprintf('%.1f:1', 10^(-Relxticks(i)));
    else
        Relxticklabels{i} = sprintf('1:%.1f', 10^(Relxticks(i)));
    end
end
Relxticklabels{1} = char(Drugs(1));
Relxticklabels{end} = char(Drugs(2));
xvals = [Relxticks(1) t_rfits.RelBAratio' Relxticks(end)];

%%%%%%%%%%%%%%%
% first plot for the CI
% conc ratio axis
get_newaxes(pos(1,:),1)
set(gca,'xtick',Relxticks,'xticklabel',xticklabels,'xticklabelrotation',90, ...
    'ytick',[],'xaxislocation','top',Plotting_parameters.axes{:})
ylim(log2([min([.45;CI(isfinite(CI))]) max([2.2;CI(isfinite(CI))])]))
xlim(Relxticks([1 end]))
xlabel('Drug ratios (Conc)',Plotting_parameters.axislabel{:})


get_newaxes(pos(1,:),1)
plot(repmat([xvals([1 end]) NaN],1,3), reshape([1;1;NaN]*[-1 0 1],[],1)', ...
    '-', 'color', [.7 .7 .7])
plot([0 0], [-32 32], '-', 'color', [.7 .7 .7])

GIcolors = [0 0 0; winter(length(extraGIs))];
h = NaN(length(GIvalues)+1,1);
for j=[2:length(GIvalues) 1] % plot GI50 last
    CIvals = [1 CI(:,j)' 1];
    h(j) = plot(xvals,log2(CIvals), 'o:','color',GIcolors(j,:),'linewidth',.5);
    h(length(GIvalues)+1) = plot(xvals(isfinite(CIvals)),log2(CIvals(isfinite(CIvals))), '--', ...
        'color',GIcolors(j,:));
    plot(allGIs(j).RelBAratio, log2(allGIs(j).CI), 'o:','color',GIcolors(j,:),...
        'linewidth',1);
end

set(gca,'xtick',Relxticks,'xticklabel',Relxticklabels,'xticklabelrotation',90, ...
    'ytick',-5:5,'yticklabel',[strcat('1/',num2cellstr(2.^(5:-1:1))) ...
    strcat(num2cellstr(2.^(0:5)),'/1')],Plotting_parameters.axes{:})
ylim(log2([min([.45;CI(isfinite(CI))]) max([2.2;CI(isfinite(CI))])]))
xlim(Relxticks([1 end]))
ylabel('Combination index',Plotting_parameters.axislabel{:})
xlabel(['Drug ratios (' strjoin(NormConc, '/') 'normalized)'],Plotting_parameters.axislabel{:})
% legend(h,[strcat('GI', num2cellstr(GIvalues)) 'Failed points'],'location','best', ...
%     Plotting_parameters.axislabel{:})


%%%%%%%%%%%%%%%%%%%
% plot for the Hill coeff
get_newaxes(pos(2,:),1)
plot([0 0], [-32 32], '-', 'color', [.7 .7 .7])
plot(xvals, [Hills(1) t_rfits.Hill' Hills(2)], '-k')
set(gca,'xtick',Relxticks,'xticklabel',Relxticklabels,'xticklabelrotation',90,...
    Plotting_parameters.axes{:})
xlim(Relxticks([1 end]))
ylabel('Hill coefficient',Plotting_parameters.axislabel{:})
xlabel('Drug ratios (IC50 normalized)',Plotting_parameters.axislabel{:})
ylim([min([.4 t_rfits.Hill']) max([4 t_rfits.Hill'])])


% plot for the Emax/inf
get_newaxes(pos(3,:),1)
plot([0 0], [-32 32], '-', 'color', [.7 .7 .7])
plot(xvals([1 end]), [0 0], '-', 'color', [.7 .7 .7])
h = plot(xvals, 1-[Einfs(1) t_rfits.Einf' Einfs(2)], '--k');
h(2) = plot(xvals, 1-[Emaxs(1) t_rfits.Emax' Emaxs(2)], '-k');
set(gca,'xtick',Relxticks,'xticklabel',Relxticklabels,'xticklabelrotation',90,...
    Plotting_parameters.axes{:})
xlim(Relxticks([1 end]))
ylim([min([-.3 1-t_rfits.Emax' 1-Emaxs]) max([.3 1-t_rfits.Emax' 1-Emaxs])])
ylabel('Lowest growth (1-E_{inf}/E_{max})',Plotting_parameters.axislabel{:})
xlabel('Drug ratios (IC50 normalized)',Plotting_parameters.axislabel{:})
legend(h,{'1-E_{inf}' '1-E_{max}'},'location','best',Plotting_parameters.axislabel{:})



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
    'xticklabel', num2cellstr(Conc{1}(2:2:end),'%.2g'),...
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

% plot the isoGI cruves
h = NaN(length(GIvalues)+1,1);
for j=[ 2:length(GIvalues) 1]
    % theoric curve (additive)
    h(j) = plot(log10(NullGIConc{2}(:,j)), log10(NullGIConc{1}(:,j)), ...
        ':','color', GIcolors(j,:), 'linewidth',1);
    % real curve
    h(j) = plot(log10(GIConc{2}(:,j)), log10(GIConc{1}(:,j)), ...
        'o:','color', GIcolors(j,:), 'linewidth',1);
    
    % theoric curve (additive)
    h(j) = plot(log10(allGIs(j).NullConc(:,2)), log10(allGIs(j).NullConc(:,1)), ...
        '-','color', GIcolors(j,:), 'linewidth',1);
    % real curve
    h(j) = plot(log10(allGIs(j).RealConc(:,2)), log10(allGIs(j).RealConc(:,1)), ...
        'o-','color', GIcolors(j,:), 'linewidth',1);
    
end
h(end) = plot([NaN NaN], [NaN NaN], '--', 'color', [.3 .3 .3]);


% plot the ratio directions
for i=1:height(t_rfits)
    t_sub = sortrows(t_combodata(t_combodata.BAratio==t_rfits.BAratio(i),:),'Conc');
    plot(log10(t_sub.Conc2), log10(t_sub.Conc),'-', 'color', [.7 .7 .7]);
end


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
