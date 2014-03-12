
% function [ICfit, GIfit] = extract_ICcurves(t_CL, Drugs, removed_replicates, savefile)

removed_replicates = {[] 3 []};
if ~exist('removed_replicates','var')
    removed_replicates = rempat({[]}, 1, length(t_CL));
else
    assert(length(removed_replicates)==length(t_CL))
end

%%
CellLines = categorical(1,length(t_CL));
for i=1:length(t_CL)
    CellLines(i) = t_CL{i}.CellLine(1);
end
%% plot the controls and the single drugs titration


allICfits = cell(1,3);

for iCL=1:length(t_CL)
    get_newfigure(iCL+10*i_corr)
    
    CtrlIdx = t_CL{iCL}.CtrlIdx==1;
    allICfits{iCL} = cell(1,length(Drugs{iCL}));
    
    GoodReplicates = ToRow(unique(t_CL{iCL}.Replicate));
    if ~isempty(removed_replicates{iCL})
        GoodReplicates = setdiff(ToRow(unique(t_CL{iCL}.Replicate)), removed_replicates{iCL});        
    end
    GoodReplicates = GoodReplicates(GoodReplicates>0);
    
    if ismember(t_CL{iCL}.Properties.VariableNames, 'Day0') && any(t_CL{iCL}.Day0)
        DoGI50 = true;
        SeededNumber = trimmean(t_CL{iCL}.Cellcount(t_CL{iCL}.Day0),50);
    else
        DoGI50 = false;
    end
    
    for iD=1:length(Drugs{iCL})
        subplot(3,3,iD)
        
        DrugCidx = find(strcmp(t_CL{iCL}.Properties.VariableNames,Drugs{iCL}(iD).DrugName));
        
        otherDrugCidx = find(ismember(t_CL{iCL}.Properties.VariableNames,...
            {Drugs{iCL}(setdiff(1:length(Drugs{iCL}),iD)).DrugName}));
        
        DrugIdx = table2array(t_CL{iCL}(:,DrugCidx))~=0 & ...
            all(table2array(t_CL{iCL}(:,otherDrugCidx))==0,2);
        Doses = setdiff(t_CL{iCL}.(Drugs{iCL}(iD).DrugName)(DrugIdx),0);
        Relcnt = NaN(max(t_CL{iCL}.Replicate), length(Doses));
        
        plot([min(Doses) max(Doses)], [1 1],'-c')
        hold on
        
        for iR = 1:max(t_CL{iCL}.Replicate)
            RepIdx = t_CL{iCL}.Replicate==iR;
            
            ctrl = mean(t_CL{iCL}.Nuclei_NumberOfObjects(CtrlIdx & RepIdx));
            std_ctrl = std(t_CL{iCL}.Nuclei_NumberOfObjects(CtrlIdx & RepIdx));
            
            if iD==1
                fprintf('%s ctrl: %.0f +/- %.0f (%i)\n', char(CellLines(iCL)), ctrl, ...
                    std_ctrl,sum(CtrlIdx & RepIdx));
            end
            
            t_doses = sortrows(t_CL{iCL}(RepIdx & DrugIdx,:),Drugs{iCL}(iD).DrugName);
            assert(all(t_doses.(Drugs{iCL}(iD).DrugName)==Doses))
            
            Relcnt(iR,:) = t_doses.Nuclei_NumberOfObjects/ctrl;
            if ismember(iR, GoodReplicates)
                plot(Doses,Relcnt(iR,:),'.-');
            else
                plot(Doses,Relcnt(iR,:),'.-','color',[.7 .7 .7]);
            end
            
            plot([min(Doses) max(Doses)], [1 1]-(std_ctrl/ctrl),':c')
            
        end
        
        plot(Doses, mean(Relcnt(GoodReplicates,:)), '.-k')
        
        [IC50, Hill, Emax, r2, fit_final, p, log] = ...
            ICcurve_fit(Doses,  mean(Relcnt(GoodReplicates,:)));
        allICfits{iCL}{iD} = fit_final;
        plot(Doses, fit_final(Doses), '.-r','linewidth',2)
        
        title(sprintf('%s, r=%.2f, %s',Drugs{iCL}(iD).DrugName, r2, tag))
        set(gca,'xscale','log')
        ylim([0 1.4])
        xlim([min([Drugs{iCL}.SingleDoses]) max([Drugs{iCL}.SingleDoses])])
    end
    
    set(gcf,'color','w','position',[50 450 900 550], ...
        'PaperUnits','centimeters','papersize',[28 18], 'PaperPositionMode', 'auto', ...
        'FileName',[folder '/ICcurves_' char(CellLines(iCL)) '_' timestamp '_' tag '.pdf'])
    
    
end

if savefile
save(['L1000validation_ICfits_' timestamp '_' tag '.mat'], 'allICfits', 'CellLines', 'Drugs', 't_CL');
end

function

