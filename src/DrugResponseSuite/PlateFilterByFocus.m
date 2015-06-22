function t_data = PlateFilterByFocus(t_data)
% t_data = PlateFilterByFocus(t_data)

t_data.filtered = false(height(t_data),1);
t_data.delta = NaN(height(t_data),1);

if ~all(isvariable(t_data, {'Column' 'Row'}))
    [Row,Column] = ConvertWellsToRowCol(cellstr(t_data.Well));
    t_data = [t_data table(Row,Column)];
end

plates = unique(t_data.Barcode);

for iP = 1:length(plates)
    
    subt = t_data(t_data.Barcode == plates(iP),:);
    
    maxrow = 2^ceil(log2(max(subt.Row)));
    maxcol = 3*2^ceil(log2(max(subt.Column)/3));
    subt.DistEdge = min([(subt.Row) (maxrow-subt.Row+1) ...
        (subt.Column) (maxcol-subt.Column+1)],[],2);
    
    [data, labels] = table_to_ndarray(subt, 'keyvars', {'Row' 'Column' 'Time'}, ...
        'outer', 1, 'valvars', {'focus' 'Cellcount'});
    
    %%
    %     get_newfigure(998,[50 550 900 400])
    
    for iT=1:length(labels{3}.Time)
    plate = data(:,:,iT,1);
    medp = median(plate(:));
        fprintf('.');
        %     filt = [1 1 1; 1 0 1; 1 1 1]/8;
        %     filt = ones(5)/25;
        %     iplate(:,:,iT) = filter2(filt, plate(:,:,iT));
        
        [x,y] = meshgrid(labels{2}.Column, labels{1}.Row);
        x = x(:); y = y(:);
        focus = reshape(plate-median(plate(:)),[],1);
        focused = (focus>quantile(focus,.25)) & (focus<quantile(focus,.75));
        plate_fit = fit([x(focused) y(focused)], focus(focused), ...
            'poly22', 'robust', 'Bisquare');
        
        
        iplate = reshape(plate_fit(x, y), size(plate));
        delta = iplate - plate;
        delta = delta - median(delta(:));
        
        diffp = diff(quantile(delta(:),[.25 .75]));
    
        %         clf
        %         subplot(131)
        %         surf(iplate(:,:,iT));
        %         hold on
        %         surf(plate(:,:,iT));
        %
        %         subplot(132)
        %         imagesc(delta, [-4 4])
        %         colormap( [1 0 0; gray(200); 0 0 1])
        %
        %         subplot(133)
        %         cnt = log10(data(:,:,iT,2));
        %
        %         plot(delta(:), cnt(:),'.')
        %         pause
        idx1 = t_data.Barcode==plates(iP) & t_data.Time==labels{3}.Time(iT);
        for iW = 1:numel(delta)
            idx = idx1 & t_data.Well==ConvertRowColToWells(x(iW), y(iW));
            t_data.filtered(idx) = (delta(iW)>2*diffp) | (delta(iW)<2*diffp);
            t_data.delta(idx) = delta(iW)*diffp;
        end
        %         if any(abs(delta(:))>1)
        %             [filt_Row, filt_Col] = find(abs(delta)>1); true;
        %         end
        
    end
    fprintf('\n');
end
