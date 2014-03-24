function h = get_subaxes(Nrows, Ncols, RowIdx, ColIdx, holded, offsets)
% h = get_subaxes(Nrows, Ncols, RowIdx, ColIdx, holded, offsets)
%   generate new axes such that it will fit (Nrows, Ncols) axes.
%
%   (Nrows, Ncols) are the number of rows and columns
%   (RowIdx, ColIdx) are the indexes for row and column
%
%   holded is optional (default not hold)
%   offsets is optional (default rather compact: [.08 .9 .06 .91])
%

if ~exist('holded','var')
    holded = false;
end
if ~exist('offsets','var')
    offsets = [.08 .9 .06 .91];
end


temph = get_newaxes([offsets(1)+offsets(2)*(ColIdx-1)/Ncols ...
    offsets(3)+offsets(4)*(RowIdx-1)/Nrows ...
    -offsets(1)+(offsets(2)-.01)/Ncols ...
    -offsets(3)+(offsets(4)-.01)/Nrows], holded);

if nargout>0
    h = temph;
end