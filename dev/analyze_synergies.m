function [a] = analyze_synergies(t_CL, t_Results, Drugs, removed_replicates, figtag)

Generate_Plotting_parameters

if ~exist('removed_replicates','var') || isempty(removed_replicates)
    removed_replicates = [];
end


assert(HasUniqueElement(t_CL.CellLine))
CellLine = char(unique(t_CL.CellLine));


%% plot the controls and the single drugs titration

removed_replicates = [];

CtrlIdx = t_CL.Ctrl==1;
meanCtrl = varfun(@mean,t_CL(CtrlIdx,:),'InputVariables', ...
    'Cellcount','GroupingVariables','Replicate');
meanCtrl = meanCtrl.mean_Cellcount;

if isfield(Drugs,'Primary')
    PrimDrugCidx = find([Drugs.IsPrimary]);
    SecDrugCidx = find(~[Drugs.IsPrimary]);
    ComboDidx = Allpairs(PrimDrugCidx, SecDrugCidx);
    Nrows = length(PrimDrugCidx);
    Ncols = length(SecDrugCidx);
else
    ComboDidx = nchoosek(1:length(Drugs),2);
    Nrows = length(Drugs)-1;
    Ncols = length(Drugs);
end

get_newfigure(figtag+2)
get_newfigure(figtag+1)

for iD = 1:size(ComboDidx,1)
    Doses1 = [0 Drugs(ComboDidx(iD,1)).ComboDoses];
    Drug1 = Drugs(ComboDidx(iD,1)).DrugName;
    Drug1CIdx = find(strcmp(t_CL.Properties.VariableNames, Drug1));
    
    Drug1Idx = table2array(t_CL(:,Drug1CIdx))~=0 ;
    
    Doses2 = [0 Drugs(ComboDidx(iD,2)).ComboDoses];
    Drug2 = Drugs(ComboDidx(iD,2)).DrugName;
    
    Drug2CIdx = find(strcmp(t_CL.Properties.VariableNames, Drug2));
    
    
    ComboIdx = table2array(t_CL(:,Drug1CIdx))~=0 & ...
        table2array(t_CL(:,Drug2CIdx))~=0;
    
    Drug1fit = t_Results.ICfit{ComboDidx(iD,1)};
    Drug2fit = t_Results.ICfit{ComboDidx(iD,2)};
    
    CellCnt = NaN(length(Doses1), length(Doses2), max(t_CL.Replicate));
    
    AddRelcnt = NaN(length(Doses1), length(Doses2));
    
    for iDo1 = 1:length(Doses1)
        for iDo2 = 1:length(Doses2)
            if iDo1 == 1 && iDo2 == 1
                CellCnt(iDo1,iDo2,:) = meanCtrl;
                AddRelcnt(iDo1,iDo2) = 1;
            elseif iDo1 == 1
                CellCnt(iDo1,iDo2,:) = meanCtrl*Drug2fit(Doses2(iDo2));
                AddRelcnt(iDo1,iDo2) = min(1,Drug2fit(Doses2(iDo2)));
            elseif iDo2 == 1
                CellCnt(iDo1,iDo2,:) = meanCtrl*Drug1fit(Doses1(iDo1));
                AddRelcnt(iDo1,iDo2) = min(1,Drug1fit(Doses1(iDo1)));
            else
                AddRelcnt(iDo1,iDo2) = ...
                    min(1,Drug1fit(Doses1(iDo1)))*...
                    min(1,Drug2fit(Doses2(iDo2)));
                for iR = setdiff(1:max(t_CL.Replicate),removed_replicates)
                    CellCnt(iDo1,iDo2,iR) = t_CL.Cellcount(...
                        t_CL.(Drug1)==Doses1(iDo1) & ...
                        t_CL.(Drug2)==Doses2(iDo2) & ...
                        t_CL.Replicate==iR);
                end
            end
        end
    end
    
    Relcnt = nanmean(CellCnt./repmat(reshape(meanCtrl,1,1,[]),...
        length(Doses1),length(Doses2)),3);
    
    for iF=1:2
        figure(figtag+iF)
        get_newaxes([.04+.92*(ComboDidx(iD,2)-1)/Ncols ...
            .05+.9*(ComboDidx(iD,1)-1)/Nrows ...
            -.06+.92/Ncols -.08+.9/Nrows])
        
        if iF==1
            imagesc(Relcnt,[.2 1])
            colormap(Plotting_parameters.cmapWP)
        else
            delta = (min(Relcnt,1)-AddRelcnt);
            Reldelta = delta./AddRelcnt;
            Reldelta(abs(delta)<.05) = 0;
            delta(abs(delta)<.05) = 0;
            
            imagesc(delta,[-.301 .301])
            
            colormap(Plotting_parameters.cmapBR)
            
        end
        if iD==1            
            hc = colorbar('location','northoutside');
            set(hc,'position',[.4 .9 .2 .02],'fontsize',8)
            if iF==1
                xlabel(hc,'Cell count (rel to ctrl)','fontsize',10,'fontweight','bold')
            elseif iF==2
                xlabel(hc,'Delta Rel cell count with prediction','fontsize',10,'fontweight','bold')
            end
        end
        
        set(gca,'xtick',1:length(Doses2),'xticklabel',1e-3*round(1e3*[0 Doses2(2:end)]), ...
            'ytick',1:length(Doses1),'yticklabel',1e-3*round(1e3*[0 Doses1(2:end)]),'fontsize',8)
        ylabel(Drug1)
        xlabel(Drug2)
    end
    
end

for iF=1:2
    figure(figtag+iF)
    set(gcf,'color','w','position',[400 30+(iF-1)*250 900 550], ...
        'PaperUnits','centimeters','papersize',[28 18], 'PaperPositionMode', 'auto')
    if iF==1
        set(gcf,'FileName',['./RelCellCnt_' CellLine '.pdf'])
    elseif iF==2
        set(gcf,'FileName',['./DeltaCombo_' CellLine '.pdf'])
    end
end




