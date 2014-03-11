timestamp = '20140220';
tag = 'corrected';

load(['L1000validation_ICfits_' timestamp '_' tag '.mat']);



CellLines = categorical(1,length(t_CL));
for i=1:length(t_CL)
    CellLines(i) = t_CL{i}.CellLine(1);
end

%% plot the controls and the single drugs titration

removed_replicates = {[] 3 []};

for iCL=1:3
    CtrlIdx = all(table2array(t_CL{iCL}(:,7:end))==0,2);
    meanCtrl = varfun(@mean,t_CL{iCL}(CtrlIdx,:),'InputVariables', ...
        'Nuclei_NumberOfObjects','GroupingVariables','Replicate');
    meanCtrl = meanCtrl.mean_Nuclei_NumberOfObjects;
    
    PrimDrugCidx = find([Drugs{iCL}.IsPrimary]);
    SecDrugCidx = find(~[Drugs{iCL}.IsPrimary]);
    
    get_newfigure(310+iCL)
    get_newfigure(320+iCL)
    get_newfigure(330+iCL)
    
    for iD1=1:length(PrimDrugCidx)
        PrimDrug = Drugs{iCL}(PrimDrugCidx(iD1)).DrugName;
        Doses1 = [1e-4 Drugs{iCL}(PrimDrugCidx(iD1)).ComboDoses];
        
        Drug1Cidx = find(strcmp(t_CL{iCL}.Properties.VariableNames,PrimDrug));
        
        Drug1Idx = table2array(t_CL{iCL}(:,Drug1Cidx))~=0 ;
        
        for iD2=1:length(SecDrugCidx)
            SecDrug = Drugs{iCL}(SecDrugCidx(iD2)).DrugName;
            Doses2 = [1e-4 Drugs{iCL}(SecDrugCidx(iD2)).ComboDoses];
            
            Drug2Cidx = find(strcmp(t_CL{iCL}.Properties.VariableNames, SecDrug));
            Drug2Idx = table2array(t_CL{iCL}(:,Drug2Cidx))~=0 ;
            
            
            ComboIdx = table2array(t_CL{iCL}(:,Drug1Cidx))~=0 & ...
                table2array(t_CL{iCL}(:,Drug2Cidx))~=0;
            
            
%             Doses1 = [1e-4 unique(table2array(t_CL{iCL}(ComboIdx,PrimDrug)))'];
%             Doses2 = [1e-4 unique(table2array(t_CL{iCL}(ComboIdx,SecDrug)))'];
            
            CellCnt = NaN(length(Doses1), length(Doses2), max(t_CL{iCL}.Replicate));
            
            AddRelcnt = NaN(length(Doses1), length(Doses2));
            
            for iDo1 = 1:length(Doses1)
                for iDo2 = 1:length(Doses2)
                    if iDo1 == 1 && iDo2 == 1
                        CellCnt(iDo1,iDo2,:) = meanCtrl;
                        AddRelcnt(iDo1,iDo2) = 1;
                    elseif iDo1 == 1
                        CellCnt(iDo1,iDo2,:) = meanCtrl*allICfits{iCL}{SecDrugCidx(iD2)}(Doses2(iDo2));
                        AddRelcnt(iDo1,iDo2) = min(1,allICfits{iCL}{SecDrugCidx(iD2)}(Doses2(iDo2)));
                    elseif iDo2 == 1
                        CellCnt(iDo1,iDo2,:) = meanCtrl*allICfits{iCL}{PrimDrugCidx(iD1)}(Doses1(iDo1));
                        AddRelcnt(iDo1,iDo2) = min(1,allICfits{iCL}{PrimDrugCidx(iD1)}(Doses1(iDo1)));
                    else
                        AddRelcnt(iDo1,iDo2) = min(1,allICfits{iCL}{PrimDrugCidx(iD1)}(Doses1(iDo1)))*...
                            min(1,allICfits{iCL}{SecDrugCidx(iD2)}(Doses2(iDo2)));
                        for iR = setdiff(1:max(t_CL{iCL}.Replicate),removed_replicates{iCL})
                            CellCnt(iDo1,iDo2,iR) = t_CL{iCL}.Nuclei_NumberOfObjects(...
                                t_CL{iCL}.(PrimDrug)==Doses1(iDo1) & ...
                                t_CL{iCL}.(SecDrug)==Doses2(iDo2) & ...
                                t_CL{iCL}.Replicate==iR);
                        end
                    end
                end
            end
            
            Relcnt = nanmean(CellCnt./repmat(reshape(meanCtrl,1,1,[]),...
                length(Doses1),length(Doses2)),3);
            
            for iF=1:3
                figure(300+10*iF+iCL)
                get_newaxes([(iD2-1)*.16+.06 (3-iD1)*.29+.07 .1 .2])
                
                if iF==1
                    imagesc(Relcnt,[.2 1])
                    colormap(Plotting_parameters.cmapWP)
                else
                    delta = (min(Relcnt,1)-AddRelcnt);
                    Reldelta = delta./AddRelcnt;
                    Reldelta(abs(delta)<.05) = 0;
                    delta(abs(delta)<.05) = 0;
                    
                    if iF==2
                        imagesc(delta,[-.301 .301])
                    else
                        imagesc(Reldelta,[-.501 .501])
                    end
                    
                    colormap(Plotting_parameters.cmapBR)
                    
                end
                if iD1==1 && iD2==1
                    hc = colorbar('location','northoutside');
                    set(hc,'position',[.4 .9 .2 .02],'fontsize',8)
                    if iF==1
                        xlabel(hc,'Cell count (rel to ctrl)','fontsize',10,'fontweight','bold')
                    elseif iF==2
                        xlabel(hc,'Delta Rel cell count with prediction','fontsize',10,'fontweight','bold')
                    else
                        xlabel(hc,'Relative Delta with prediction','fontsize',10,'fontweight','bold')
                    end
                end
                
                set(gca,'xtick',1:length(Doses2),'xticklabel',[0 Doses2(2:end)], ...
                    'ytick',1:length(Doses1),'yticklabel',[0 Doses1(2:end)],'fontsize',8)
                ylabel(PrimDrug)
                xlabel(SecDrug)
            end
        end
    end
    
    for iF=1:3
        figure(300+10*iF+iCL)
        set(gcf,'color','w','position',[-350+iCL*400 30+(iF-1)*250 900 550], ...
            'PaperUnits','centimeters','papersize',[28 18], 'PaperPositionMode', 'auto')
        if iF==1
            set(gcf,'FileName',['./PlotsCombo/RelCellCnt_' char(CellLines(iCL)) '.pdf'])
        elseif iF==2
            set(gcf,'FileName',['./PlotsCombo/DeltaCombo_' char(CellLines(iCL)) '.pdf'])
        elseif iF==3
            set(gcf,'FileName',['./PlotsCombo/RelDeltaCombo_' char(CellLines(iCL)) '.pdf'])
        end
    end
end


