
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
t_normdata = t_data(:,{'DrugName' 'Conc' 'DrugName2' 'Conc2' 'RelGrowth'});
t_normdata.RelConc = NaN(height(t_normdata),1);
t_normdata.RelConc2 = NaN(height(t_normdata),1);
t_normdata.RelConc2(t_normdata.Conc2==0) = 0;
for i=1:2
    t_temp = sortrows(t_data(t_data.DrugName==Drugs(i) & t_data.DrugName2=='-',:), ...
        'Conc');
    Emaxs(i) = 1-min(t_temp.RelGrowth);
    
    opt = struct('pcutoff',.1,'extrapolrange',30, 'plotting',1);
    [GIs(1,i), Hills(i), Einfs(i), ~, r2s(i), EC50s(i), fitfct] = ...
        ICcurve_fit(t_temp.Conc, t_temp.RelGrowth, 'GI50', opt);
    %%%% replace by the actual subfunction and stack it in the file
    %%%%%%%%%%
    
    fprintf('\tfit for %-12s: EC50=%.3g, GI50=%.3g, r2=%.2f\n', char(Drugs(i)), EC50s(i), GIs(1,i), r2s(i));
    warnassert(r2s(i)>.7, 'Bad fit !')
    %%%%% should that stop the algorithm ??
    
    % evaluation of the other GIs
    fitopts = optimoptions('fsolve', 'Display','none');
    for j=1:length(extraGIs)
        if Einfs(i)>extraGIs(j)/100
            GIs(1+j,i) = fsolve(@(x) fitfct(x)-(1-extraGIs(j)/100), ...
                median(t_temp.Conc), fitopts);
        end
    end
    %%%%%% the other GI should directly be evaluated in the ICcurve_fit
    %%%%%% function and capped as the GI50 if the values are extrapolated
    
    %%%%%%%%%%%%%%
    %%%%  how to deal with normalization when there is no GI50 for normalization ????
    %%%%  use EC50??? instead (with a warning?)
    %%%%%%%%%%%%%%
    if isfinite(GIs(1,i))
        NormConc{i} = 'GI50';
        NormFactor(i) = GIs(1,i);
    else
        NormConc{i} = 'EC50';
        NormFactor(i) = EC50s(i);
        warnprintf('Normalization by EC50 instead of GI50')
    end
    
    idx = t_data.DrugName==Drugs(i);
    t_normdata.RelConc(idx) = t_normdata.Conc(idx)/NormFactor(i);
    idx = t_data.DrugName2==Drugs(i);
    t_normdata.RelConc2(idx) = t_normdata.Conc2(idx)/NormFactor(i);
end

assert(all(~isnan(t_normdata.RelConc)))
assert(all(~isnan(t_normdata.RelConc2)))

%% caculate the ratios on which to fit the sigmoidal cruves

t_temp = t_normdata(t_normdata.DrugName==Drugs(1) & t_normdata.DrugName2==Drugs(2),...
    {'Conc' 'Conc2' 'RelConc' 'RelConc2' 'RelGrowth'});

BAratio = 1e-2*round(log10(t_temp.RelConc2./t_temp.RelConc)*1e2);
t_temp = [t_temp table(BAratio) ...
    table(t_temp.RelConc+t_temp.RelConc2, 'variablename', {'RelConcTot'}) ];
BAratio = unique(t_temp.BAratio,'sorted');
ratiocounts = hist(t_temp.BAratio, BAratio);

BAratio = BAratio(ratiocounts>4);
fprintf('\tFitting for %i ratios:\n', length(BAratio));

t_rfits = [table(BAratio) array2table(NaN(length(BAratio),5+length(extraGIs)), ...
    'variablenames', ...
    [{'GI50' 'Hill' 'Einf' 'r2' 'Emax'} strcat('GI', num2cellstr(extraGIs))])];
colors = parula(length(BAratio));

figure(999);clf;hold on;
for i = 1:height(t_rfits)
    t_sub = sortrows(t_temp(t_temp.BAratio == t_rfits.BAratio(i),:), ...
        'RelConcTot');
    t_rfits.Emax(i) = 1-min(t_sub.RelGrowth);
    
    %%%%%%%
    opt = struct('pcutoff',1,'extrapolrange',30);
    [t_rfits.GI50(i), t_rfits.Hill(i), t_rfits.Einf(i), ~, t_rfits.r2(i), ...
        ~, fitfct, ~, log] = ...
        ICcurve_fit(t_sub.RelConcTot, t_sub.RelGrowth, 'GI50', opt);
    %%%% replace by the actual subfunction and stack it in the file
    %%%%%%%%%%
    
    for j=extraGIs
        fitopts = optimoptions('fsolve', 'Display','none');
        temp = fsolve(@(x) fitfct(x)-(1-j/100), ...
            median(t_sub.RelConcTot), fitopts);
        % capping the values to avoid too much extrapolation
        if temp<1e-2*min(t_sub.RelConcTot)
            t_rfits.(['GI' num2str(j)])(i) = -Inf;
        elseif temp>1e2*max(t_sub.RelConcTot)
            t_rfits.(['GI' num2str(j)])(i) = Inf;
        else
            t_rfits.(['GI' num2str(j)])(i) = temp;
        end
    end
    %%%%%% the other GI should directly be evaluated in the ICcurve_fit
    %%%%%% function and capped as the GI50 if the values are extrapolated
    
    fprintf('\t\tfit for BAratio %-5g: GI50=%.3g, r2=%.2f\n', t_rfits.BAratio(i), ...
        t_rfits.GI50(i), t_rfits.r2(i));
    
    plot(t_sub.RelConcTot, t_sub.RelGrowth, 'o', 'color', colors(i,:))
    
    if isinf(t_rfits.GI50(i))
        fprintf(['\t\t --> ' log '\n'])
    else
        x = 10.^(log10(min([t_rfits.GI50(i);t_sub.RelConcTot])):.01:...
            log10(max([t_rfits.GI50(i);t_sub.RelConcTot])));
        plot(x, fitfct(x), '-', 'color', colors(i,:))
        plot(x(end)*[1 1.2], [1 1]-t_rfits.Emax(i), ':', 'color', colors(i,:))
        plot(x(end)*[1 1.2], [1 1]-t_rfits.Einf(i), '--', 'color', colors(i,:))
    end
end
set(gca,'xscale','log')


%%%%%  correction for GIxx as a linear interpolation; 
w = [10.^-t_rfits.BAratio 10.^t_rfits.BAratio];
for j=extraGIs
    relGI = GIs(GIvalues==j,:)./GIs(1,:);    
    t_rfits.(['GI' num2str(j)]) = t_rfits.(['GI' num2str(j)]) ./ ...
        (w*relGI'./sum(w,2));
end

% getting the actual Concentration for the GIxx for plotting later
w = [10.^-t_rfits.BAratio 10.^t_rfits.BAratio];
for j=GIvalues
    % additive case
    t_rfits.(['Conc' num2str(j) '_1ctrl']) = w(:,1)*GIs(GIvalues==j,1)./sum(w,2);
    t_rfits.(['Conc' num2str(j) '_2ctrl']) = w(:,2)*GIs(GIvalues==j,2)./sum(w,2);
    
    % actual data
    t_rfits.(['Conc' num2str(j) '_1']) =  t_rfits.(['GI' num2str(j)]) .* ...
        w(:,1)*GIs(GIvalues==j,1)./sum(w,2);
    t_rfits.(['Conc' num2str(j) '_2']) = t_rfits.(['GI' num2str(j)]) .* ...
        w(:,2)*GIs(GIvalues==j,2)./sum(w,2);    
end

%%


get_newfigure(101,[50 50 700 550], [char(unique(t_data.CellLine)) ...
    '_Combo_' char(Drugs(1)) '+' char(Drugs(2)) '.pdf'])

% spacing for checkboard
xwidth = .03;
ywidth = .04;
space = .005;

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


% define the ratio in the normalized space
xticks = [2*t_rfits.BAratio(1)-t_rfits.BAratio(2)
    t_rfits.BAratio(2:floor(height(t_rfits)/5):(end-1))
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
xvals = [xticks(1) t_rfits.BAratio' xticks(end)];

get_newaxes(pos(1,:),1)
plot(repmat([xvals([1 end]) NaN],1,3), reshape([1;1;NaN]*[-1 0 1],[],1)', ...
    '-', 'color', [.7 .7 .7])
plot([0 0], [-32 32], '-', 'color', [.7 .7 .7])

GIcolors = [0 0 0; winter(length(extraGIs))];
h = NaN(length(GIvalues)+1,1);
allGIvals = [];
for j=[2:length(GIvalues) 1] % plot GI50 last
    GIvals = [1 t_rfits.(['GI' num2str(GIvalues(j))])' 1];
    h(j) = plot(xvals,log2(GIvals), '-','color',GIcolors(j,:),'linewidth',1);
    h(length(GIvalues)+1) = plot(xvals(isfinite(GIvals)),log2(GIvals(isfinite(GIvals))), '--', ...
        'color',GIcolors(j,:));
    allGIvals = [allGIvals; GIvals(isfinite(GIvals))'];
end

set(gca,'xtick',xticks,'xticklabel',xticklabels,'xticklabelrotation',90, ...
    'ytick',-5:5,'yticklabel',[strcat('1/',num2cellstr(2.^(5:-1:1))) ...
    strcat(num2cellstr(2.^(0:5)),'/1')],Plotting_parameters.axes{:})
ylim(log2([min([.45;allGIvals]) max([2.2;allGIvals])]))
xlim(xticks([1 end]))
ylabel('Combination index',Plotting_parameters.axislabel{:})
xlabel('Drug ratios (IC50 normalized)',Plotting_parameters.axislabel{:})
legend(h,[strcat('GI', num2cellstr(GIvalues)) 'Failed points'],'location','best', ...
    Plotting_parameters.axislabel{:})

%
get_newaxes(pos(2,:),1)
plot([0 0], [-32 32], '-', 'color', [.7 .7 .7])
plot(xvals, [Hills(1) t_rfits.Hill' Hills(2)], '-k')
set(gca,'xtick',xticks,'xticklabel',xticklabels,'xticklabelrotation',90,...
    Plotting_parameters.axes{:})
xlim(xticks([1 end]))
ylabel('Hill coefficient',Plotting_parameters.axislabel{:})
xlabel('Drug ratios (IC50 normalized)',Plotting_parameters.axislabel{:})
ylim([min([.4 t_rfits.Hill']) max([4 t_rfits.Hill'])])

get_newaxes(pos(3,:),1)
plot([0 0], [-32 32], '-', 'color', [.7 .7 .7])
plot(xvals([1 end]), [0 0], '-', 'color', [.7 .7 .7])
h = plot(xvals, 1-[Einfs(1) t_rfits.Einf' Einfs(2)], '--k');
h(2) = plot(xvals, 1-[Emaxs(1) t_rfits.Emax' Emaxs(2)], '-k');
set(gca,'xtick',xticks,'xticklabel',xticklabels,'xticklabelrotation',90,...
    Plotting_parameters.axes{:})
xlim(xticks([1 end]))
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
Conc1 = unique(t_data.Conc(t_data.DrugName==Drugs(1)));
logC1 = log10(Conc1);
Conc2 = unique(t_data.Conc(t_data.DrugName==Drugs(2)));
logC2 = log10(Conc2);
GrowthMx = NaN(length(Conc1)+1, length(Conc2)+1);
GrowthMx(1,1) = 1;

for iC1 = 1:length(Conc1)
    GrowthMx(iC1+1,1) = t_data.RelGrowth(t_data.DrugName==Drugs(1) & ...
        t_data.Conc==Conc1(iC1) & t_data.DrugName2=='-');
    
    for iC2 = 1:length(Conc2)
        if iC1==1
            GrowthMx(1,iC2+1) = t_data.RelGrowth(t_data.DrugName==Drugs(2) & ...
                t_data.Conc==Conc2(iC2) & t_data.DrugName2=='-');
        end
        
        val = t_data.RelGrowth(t_data.DrugName==Drugs(1) & ...
            t_data.Conc==Conc1(iC1) & t_data.DrugName2==Drugs(2) & ...
            t_data.Conc2==Conc2(iC2) );
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
imagesc(ones(1,length(Conc1)), logC1, GrowthMx(2:end,1),[-1 1]);
ylabel([char(Drugs(1))  '[uM]'], Plotting_parameters.axislabel{:})
set(gca,'xtick',[], 'ytick', logC1(2:2:end), ...
    'yticklabel', num2cellstr(Conc1(2:2:end),'%.2g'), ...
    Plotting_parameters.axes{:},'ydir','normal')
% plot the GIxx for Drug 1
for i=1:length(GIvalues)
    plot([.5 1.5], log10(GIs(i,1))*[1 1],'-','color', GIcolors(i,:), ...
        'linewidth',2);
end
ylim(ylims)

% plot for Drug 2
get_newaxes(subpos(2,:),1)
imagesc(logC2, ones(1,length(Conc2)), GrowthMx(1,2:end),[-1 1]);
xlabel([char(Drugs(2)) '[uM]'], Plotting_parameters.axislabel{:})
set(gca,'ytick',[], 'xtick', logC2(2:2:end), ...
    'xticklabel', num2cellstr(Conc2(2:2:end),'%.2g'),...
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
for j=1:length(GIvalues)
    % theoric curve (additive)
    plot(log10(t_rfits.(['Conc' num2str(GIvalues(j)) '_2ctrl'])), ...
        log10(t_rfits.(['Conc' num2str(GIvalues(j)) '_1ctrl'])), ...
        '--','color', GIcolors(j,:), 'linewidth',1);    
    % real curve
    h(j) = plot(log10(t_rfits.(['Conc' num2str(GIvalues(j)) '_2'])), ...
        log10(t_rfits.(['Conc' num2str(GIvalues(j)) '_1'])), ...
        '-','color', GIcolors(j,:), 'linewidth',2);
end
h(end) = plot([NaN NaN], [NaN NaN], '--', 'color', [.3 .3 .3]);
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
