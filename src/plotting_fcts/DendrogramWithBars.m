function axis_handles = DendrogramWithBars(SampleDist, t_SampleLabels, AxisPos, xwidth, ...
  PlotOptions, colorthreshold)
% function axis_handles = DendrogramWithBars(SampleDist, t_SampleLabels, AxisPos, xwidth, ...
%   PlotOptions, colorthreshold)

labels = cellstr2str( table2cellstr(t_SampleLabels,0));


axis_handles = get_newaxes(AxisPos);
tree = linkage(SampleDist,'average');
if ~exist('colorthreshold','var') || isempty(colorthreshold)
    [h,~,outperm] = dendrogram(tree,0,'labels',labels,'orientation','left');
    for i=1:length(h), set(h(i),'color','k'),end
else
    [h,~,outperm] = dendrogram(tree,0,'labels',labels,'orientation','left',...
        'colorthreshold',colorthreshold);
    SampleColor = NaN(length(labels),3);
    for i=1:length(h)
        idx = find(get(h(i),'xdata')==0);
        if ~isempty(idx)
            ypos = get(h(i),'ydata');
            for j=1:length(idx)
                if mod(ypos(idx(j)),1)==0
                    SampleColor(ypos(idx(j)),:) = get(h(i),'color');
                end
            end
        end
    end
end

ylabs = get(gca,'yticklabel');
set(gca,'fontsize',8, 'ytick',[])
ylim([.5 length(labels)+.5])
xlim([0 1.1])

%%
fields = varnames(t_SampleLabels);

for iF = 1:length(fields)
    axis_handles(iF+1) = get_newaxes([AxisPos(1)+AxisPos(3)+(iF-1)*xwidth AxisPos(2) xwidth AxisPos(4)],1);
    
    for i=1:length(PlotOptions.(fields{iF}))
        temp = t_SampleLabels.(fields{iF});
        temp = temp(outperm);
        idx = find(temp==PlotOptions.(fields{iF})(i));
        barh([-2;-1;idx], ones(size(idx,1)+2,1), 1, 'facecolor', ...
            PlotOptions.([fields{iF} 'Colors'])(i,:), ...
            'linestyle','none')
    end
    
    for i = find(any(diff(SampleColor,1)~=0,2))'
        plot([0 1], [i i]+.5, '-k')
    end
    
    ylim([.5 length(labels)+.5])
    xlim([.1 .9])
    set(gca,'xtick',[],'ytick',[])
end
set(gca,'ytick',1:length(labels), 'yticklabel', ylabs, 'fontsize',6, ...
    'yaxislocation','right','xtick',[])
