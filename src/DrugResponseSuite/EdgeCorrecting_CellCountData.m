function t_corrected = EdgeCorrecting_CellCountData(t_annotated, numericfields)
fprintf('Edge correction function:\n');

if ~exist('numericfields','var')
    labelfields = [{'pert_type' 'RelCellCnt' 'RelGrowth' 'DesignNumber' 'Barcode' ...
        'Untrt' 'Cellcount' 'Date' 'Row' 'Column' 'Well' 'TreatmentFile' 'Replicate', ...
        'CellLine' 'Time' 'DrugName' 'Conc'} ...
        strcat('Conc', cellfun(@(x) {num2str(x)}, num2cell(2:9)))];
    numericfields = setdiff(t_annotated.Properties.VariableNames( ...
        all(cellfun(@isnumeric, table2cell(t_annotated)))), labelfields);
    if ~isempty(numericfields)
        fprintf('\tThese numeric fields will be corrected (use ''numericfields'' to specify):\n');
        for i=1:length(numericfields)
            fprintf('\t - %s\n', numericfields{i});
        end
    end
end


barcodes = unique(t_annotated.Barcode);
t_corrected = t_annotated;

for ib = 1:length(barcodes)
    idx = find(t_annotated.Barcode==barcodes(ib));
    t_plateinfo = t_annotated(idx(1),{'Barcode','CellLine','Time'});
    t_corrected(idx,:) = EdgeCorrection(t_annotated(idx,:), ...
        numericfields, t_plateinfo);
end

end


function t_out = EdgeCorrection(t_in, numericfields, t_plateinfo)

t_out = t_in;

if ~all(isvariable(t_in, {'Col' 'Row'}))
    [Row,Col] = ConvertWellsToRowCol(cellstr(t_in.Well));
    t_in = [t_in table(Row,Col)];
end

t_ctrl = t_in(t_in.pert_type=='ctl_vehicle',:);

if any(t_ctrl.Row==1) && any(t_ctrl.Col==1) && ...
        mod(log2(max(t_ctrl.Row)),1)==0 && mod(log2(max(t_ctrl.Col/3)),1)==0
    
    edge_ctrl_idx = t_ctrl.Row==1 | t_ctrl.Col==1 | t_ctrl.Row==max(t_in.Row) ...
        | t_ctrl.Col==max(t_in.Col);
    edge_idx = t_in.Row==1 | t_in.Col==1 | t_in.Row==max(t_in.Row) ...
        | t_in.Col==max(t_in.Col);
    
    for iFields = [{'Cellcount'} numericfields]
        edge_vals = t_ctrl.(iFields{:})(edge_ctrl_idx);
        center_vals = t_ctrl.(iFields{:})(~edge_ctrl_idx);
        ratio = mean(center_vals)/mean(edge_vals);
        pval = ranksum(edge_vals, center_vals);
        % correct for edges only in their is a significant difference in
        % the controls
        if pval<0.05
            fprintf('\tCorrecting the edge effect for %s in %s (ratio=%.2f, p=%.3f)\n', ...
                iFields{:}, strjoin(table2cellstr(t_plateinfo,0),'/'), ratio, pval )
            t_out.(iFields{:})(edge_idx) = t_in.(iFields{:})(edge_idx)*ratio;
        else
            fprintf('\tNo edge correctionfor %s in %s (ratio=%.2f, p=%.3f)\n', ...
                iFields{:}, strjoin(table2cellstr(t_plateinfo,0),'/'), ratio, pval )
        end
    end
end


end


