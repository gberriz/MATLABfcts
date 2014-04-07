% t_Results = extract_ICcurves(t_CL, Drugs, removed_replicates, plotting)
%
%
%


function t_Results = extract_ICcurves(t_CL, Drugs, removed_replicates, fignum, figtag)


if ~exist('removed_replicates','var') || isempty(removed_replicates)
    removed_replicates = [];
end


if ~exist('fignum','var') || isempty(fignum)
    fignum = 0;
end

if ~exist('figtag','var') || isempty(figtag)
    figtag = '';
end

%%
assert(HasUniqueElement(t_CL.CellLine))
CellLine = char(unique(t_CL.CellLine));

%% plot the controls and the single drugs titration


CtrlIdx = t_CL.Ctrl==1;

GoodReplicates = ToRow(unique(t_CL.Replicate));
if ~isempty(removed_replicates)
    GoodReplicates = setdiff(ToRow(unique(t_CL.Replicate)), removed_replicates);
end
GoodReplicates = GoodReplicates(GoodReplicates>0);

if ismember('Day0', t_CL.Properties.VariableNames) && any(t_CL.Day0)
    DoGI50 = true;
    SeededNumber = trimmean(t_CL.Cellcount(t_CL.Day0==1),50);
else
    DoGI50 = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     remove the edges if there are made of controls and there are 
%%%         enough controls in the center
%%%%%%%%%%%%%%%%

if fignum, get_newfigure(fignum), end
fprintf('Fit for %s:\n',CellLine);
ctrls = NaN(1,max(t_CL.Replicate));
for iR = 1:max(t_CL.Replicate)
    RepIdx = t_CL.Replicate==iR;
    
    ctrls(iR) = mean(t_CL.Cellcount(CtrlIdx & RepIdx));
    std_ctrl = std(t_CL.Cellcount(CtrlIdx & RepIdx));
    
    fprintf('\tUntreated ctrl (rep %i): %.0f +/- %.0f (%i)\n', ...
        iR, ctrls(iR), std_ctrl,sum(CtrlIdx & RepIdx));
    
end



if DoGI50
    fprintf('\tSeeding number: %.0f +/- %.0f\n', SeededNumber, ...
        std(t_CL.Cellcount(t_CL.Day0==1)));
    GI50 = NaN(1, length(Drugs));
    if fignum, get_newfigure(fignum+1), end
    t_GI = table;
end

t_IC = table;
t_labels = table;

for iD=1:length(Drugs)
    DrugName = Drugs(iD).DrugName;
    t_labels = [t_labels;
        table({CellLine}, {DrugName}, 'VariableNames', {'CellLine' 'Drug'})];
    
    DrugCidx = strcmp(t_CL.Properties.VariableNames,DrugName);
    otherDrugCidx = ismember(t_CL.Properties.VariableNames,...
        {Drugs(setdiff(1:length(Drugs),iD)).DrugName});
    
    DrugIdx = table2array(t_CL(:,DrugCidx))~=0 & ...
        all(table2array(t_CL(:,otherDrugCidx))==0,2);
    
    Doses = unique(t_CL.(DrugName)(t_CL.Replicate==1 & DrugIdx));
    assert(all(Doses==unique(setdiff(t_CL.(DrugName)(DrugIdx),0))))
    
    Relcnt = NaN(max(t_CL.Replicate), length(Doses));
    
    if DoGI50
        Relgrowth = NaN(max(t_CL.Replicate), length(Doses));
    end
    if fignum
        figure(fignum)
        subplot(floor(sqrt(length(Drugs))), ...
            ceil(length(Drugs)/floor(sqrt(length(Drugs)))),iD)
        plot([min(Doses) max(Doses)], [1 1],'-c')
        hold on
        
        if DoGI50
            figure(fignum+1)
            subplot(floor(sqrt(length(Drugs))), ...
                ceil(length(Drugs)/floor(sqrt(length(Drugs)))),iD)
            plot([min(Doses) max(Doses)], [1 1],'-c')
            hold on
        end
    end
    
    for iR = 1:max(t_CL.Replicate)
        RepIdx = t_CL.Replicate==iR;
        
        ctrl = mean(t_CL.Cellcount(CtrlIdx & RepIdx));
        std_ctrl = std(t_CL.Cellcount(CtrlIdx & RepIdx));
        
        t_doses = sortrows(t_CL(RepIdx & DrugIdx,:),DrugName);
        assert(all(unique(t_doses.(DrugName))==Doses))
        
        Relcnt(iR,:) = varfun(@mean,t_doses,'GroupingVar',DrugName, ...
            'InputVar', 'Cellcount', 'outputformat','uniform')/ctrls(iR);
        if fignum
            figure(fignum)
            if ismember(iR, GoodReplicates)
                plot(Doses,Relcnt(iR,:),'.-');
            else
                plot(Doses,Relcnt(iR,:),'.-','color',[.7 .7 .7]);
            end
            
            plot([min(Doses) max(Doses)], [1 1]-(std_ctrl/ctrl),':c')
            plot([min(Doses) max(Doses)], [1 1]+(std_ctrl/ctrl),':c')
        end
        
        if DoGI50
            Relgrowth(iR,:) = (varfun(@mean,t_doses,'GroupingVar',DrugName,...
                'InputVar', 'Cellcount', 'outputformat','uniform')-SeededNumber)/...
                (ctrls(iR)-SeededNumber);
            
            if fignum
                figure(fignum+1)
                if ismember(iR, GoodReplicates)
                    plot(Doses,Relgrowth(iR,:),'.-');
                else
                    plot(Doses,Relgrowth(iR,:),'.-','color',[.7 .7 .7]);
                end
                
                plot([min(Doses) max(Doses)], [1 1]-(std_ctrl/ctrl),':c')
                plot([min(Doses) max(Doses)], [1 1]+(std_ctrl/ctrl),':c')
            end
        end
    end
    
    
    [IC50, Hill, Emax, Area, r2, fit, p, log] = ...
        ICcurve_fit(Doses, mean(Relcnt(GoodReplicates,:)), 'IC50');
    
    t_IC = [t_IC;
        cell2table({IC50, Hill, Emax, Area, r2, fit, log, {Doses'}, ...
        {mean(Relcnt(GoodReplicates,:))}}, 'VariableNames', ...
        {'IC50' 'Hill' 'Emax' 'Area' 'r2' 'ICfit' 'log' 'Doses' 'Rel_CellCnt'})];
    
    if fignum
        figure(fignum)
        plot(Doses, mean(Relcnt(GoodReplicates,:)), '.-k')
        plot(Doses, fit(Doses), '.-r','linewidth',2)
        
        title(sprintf('%s - %s, r=%.2f',CellLine, DrugName, r2))
        set(gca,'xscale','log')
        ylim([0 1.4])
        xlim([min([Drugs.SingleDoses]) max([Drugs.SingleDoses])])
        if DoGI50
            plot([min(Doses) max(Doses)], [1 1]*SeededNumber/ctrl,':k')
        end
    end
    
    if DoGI50
        [GI50, Hill, Emax, Area, r2, fit, p, log] = ...
            ICcurve_fit(Doses, mean(Relgrowth(GoodReplicates,:)), 'GI50');
        
        
        t_GI = [t_GI;
            cell2table({GI50, Hill, Emax, Area, r2, fit, log, ...
            {mean(Relgrowth(GoodReplicates,:))} SeededNumber*ones(size(r2))},...
            'VariableNames', ...
            {'GI50' 'Hill_GI' 'Emax_GI' 'Area_GI' 'r2_GI' 'GIfit' ...
            'log_GI' 'Rel_growth' 'SeededNumber'})];
        
        if fignum
            figure(fignum+1)
            plot(Doses, mean(Relgrowth(GoodReplicates,:)), '.-k')
            plot(Doses, fit(Doses), '.-r','linewidth',2)
            
            
            title(sprintf('%s - %s, r=%.2f',CellLine, DrugName, r2))
            set(gca,'xscale','log')
            ylim([min([-.1 min(Relgrowth(GoodReplicates,:))]) 1.4])
            xlim([min([Drugs.SingleDoses]) max([Drugs.SingleDoses])])
        end
    end
end


if DoGI50
    t_Results = [t_labels t_IC t_GI];
else
    t_Results = [t_labels t_IC];
end

if fignum
    figure(fignum)
    set(gcf,'color','w','position',[50 450 900 550], ...
        'PaperUnits','centimeters','papersize',[28 18], 'PaperPositionMode', 'auto', ...
        'FileName',['ICcurves_' CellLine figtag '.pdf'])
    
    if DoGI50
        figure(fignum+1)
        set(gcf,'color','w','position',[50 450 900 550], ...
            'PaperUnits','centimeters','papersize',[28 18], 'PaperPositionMode', 'auto', ...
            'FileName',['GIcurves_' CellLine figtag '.pdf'])
    end
    
    
end

