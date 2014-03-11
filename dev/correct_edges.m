function t_CLcorrected = correct_edges(t_CL, CellLines)


t_CLcorrected = t_CL;

for iCL=1:length(t_CL)
    
    bins = [min(t_CL{iCL}.Cellcount) max(t_CL{iCL}.Cellcount)];
    bins = [floor(.95*bins(1)):floor(diff(bins)/30):ceil(1.05*bins(2))];
    
    
    get_newfigure(200+iCL,[100 100 700 450])
    for iR = 1:max(t_CL{iCL}.Replicate)
        t_ctrl_plate = t_CL{iCL}(t_CL{iCL}.Replicate==iR & t_CL{iCL}.Ctrl==1,...
            {'Row','Column','Cellcount'});
        
        a1 = get_newaxes([.07+(iR-1)*.33 .4 .26 .5],1);
        a2 = get_newaxes([.07+(iR-1)*.33 .05 .26 .3],1);
        
        compare_edge_center(t_ctrl_plate, bins, a1, a2, CellLines{iCL})
    end
    
    
    get_newfigure(300+iCL,[100 500 700 450])
    for iD0 = 1:max(t_CL{iCL}.Day0)
        t_ctrl_plate = t_CL{iCL}(t_CL{iCL}.Day0==iD0,{'Row','Column','Cellcount'});
        
        a1 = get_newaxes([.07+(iD0-1)*.33 .4 .26 .5],1);
        a2 = get_newaxes([.07+(iD0-1)*.33 .05 .26 .3],1);
        
        compare_edge_center(t_ctrl_plate, bins, a1, a2, CellLines{iCL})
    end
    
    get_newfigure(400+iCL,[100 500 700 450])
    for iUnt = 1:max(t_CL{iCL}.Untrt)
        t_ctrl_plate = t_CL{iCL}(t_CL{iCL}.Untrt==iUnt,{'Row','Column','Cellcount'});
        
        a1 = get_newaxes([.07+(iUnt-1)*.33 .4 .26 .5],1);
        a2 = get_newaxes([.07+(iUnt-1)*.33 .05 .26 .3],1);
        
        compare_edge_center(t_ctrl_plate, bins, a1, a2, CellLines{iCL})
    end
end

end

function compare_edge_center(t_ctrl_plate, bins, a1, a2, CellLine)

Generate_Plotting_parameters

edgeCtrl = t_ctrl_plate( t_ctrl_plate.Row==1 | t_ctrl_plate.Row==16 | ...
    t_ctrl_plate.Column==1 | t_ctrl_plate.Column==24, :);
if isempty(edgeCtrl)
    return
end

centerCtrl = setdiff(t_ctrl_plate,edgeCtrl);

set(gcf,'CurrentAxes',a1)

n = ksdensity(edgeCtrl.Cellcount,bins,'width',diff(bins([1 2]))*3);
h = plot(bins,n);

n = ksdensity(centerCtrl.Cellcount,bins,'width',diff(bins([1 2]))*3);
h(2) = plot(bins,n,'r');

legend(h, {['Edge, n=' num2str(height(edgeCtrl))] ...
    ['Center, n=' num2str(height(centerCtrl))]})
xlim([min(bins) max(bins)])
set(gca,'ytick',[])
pval = ranksum(centerCtrl.Cellcount,edgeCtrl.Cellcount);
ratio =  median(centerCtrl.Cellcount)/median(edgeCtrl.Cellcount);
title({ [CellLine '; CV=' ...
    num2str(std(t_ctrl_plate.Cellcount)/mean(t_ctrl_plate.Cellcount),'%.2f')];
    sprintf('\\Delta=%.1f, ratio=%.3f, p=%.3f', ...
    median(centerCtrl.Cellcount)-median(edgeCtrl.Cellcount), ...
    ratio, pval)});

% edgeIdx = (t_CL{iCL}.Row==1 | t_CL{iCL}.Row==16 | ...
%     t_CL{iCL}.Column==1 | t_CL{iCL}.Column==24) & t_CL{iCL}.Replicate==iR;
% 
% t_CLcorrected{iCL}.Cellcount(edgeIdx) = ...
%     t_CLcorrected{iCL}.Cellcount(edgeIdx)*ratio;

set(gca','fontsize',8)


[Ctrl2D, labels] = table_to_ndarray(t_ctrl_plate,...
    'KeyVars', {'Row' 'Column'}, 'ValVars', {'Cellcount'}, 'outer',1);
order = sortidx(labels{2}.Column);

set(gcf,'CurrentAxes',a2)

imagesc(Ctrl2D(:,order), [min(bins) max(bins)])
xlim([.5 height(labels{2})+.5])
ylim([.5 height(labels{1})+.5])

colormap(Plotting_parameters.cmapWP)
end
