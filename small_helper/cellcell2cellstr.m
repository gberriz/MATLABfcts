function cellarray = cellcell2cellstr(cellcells, idx)
%
%   cellarray = cellcell2cellstr(cellcells, idx)
%
%       convert the form  { {x}, {y}, '', {z1 z2} } to {x y '' z1z2}
%       usefully for the ouptput of regexp when applied on cell array
%

cellarray = repmat({''}, size(cellcells));
notempty = ~cellfun(@isempty, cellcells);
notempty(cellfun(@iscell,cellcells)) = ...
    ~cellfun(@(x) isequal(x, {''}), cellcells(cellfun(@iscell,cellcells)));
if ~exist('idx','var')
    cellarray(notempty) = cellfun2(@(x) [x{:}], cellcells(notempty));
else
    cellarray(notempty) = cellfun2(@(x) [x{idx}], cellcells(notempty));
end

cellarray(cellfun(@iscell,cellarray)) = cellfun2(@(x) [x{:}], cellarray(cellfun(@iscell,cellarray)));
cellarray(cellfun(@isempty,cellarray)) = {''};
