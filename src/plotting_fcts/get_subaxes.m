function h = get_subaxes(nRows, nCols, RowIdx, ColIdx, holded, varargin)
% h = get_subaxes(Nrows, Ncols, RowIdx, ColIdx, holded, varargin)
%   generate new axes such that it will fit (Nrows, Ncols) axes.
%
%   (Nrows, Ncols) are the number of rows and columns
%   (RowIdx, ColIdx) are the indexes for row and column
%
%   holded is optional (default not hold)
%   optional fields are:
%       - xspacing
%       - yspacing
%       - any name/parameter pairs valid as axis properties
%

if ~exist('holded','var') || isempty(holded)
    holded = false;
end

p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'xspacing', .05, @isnumeric)
addParameter(p, 'yspacing', .08, @isnumeric)

parse(p,varargin{:})
extra = p.Unmatched;
p = p.Results;

extravars = {};
extrafields = fieldnames(extra);
for i=1:length(extrafields)
    extravars{2*i-1} = extrafields{i};
    extravars{2*i} = extra.(extrafields{i});
end

xspacing = p.xspacing;
axis_width = (1-(nCols+1)*xspacing)/nCols;
yspacing = p.yspacing;
axis_height = (1-(nRows+1)*yspacing)/nRows;

temph = get_newaxes([xspacing*1.5+(ColIdx-1)*(axis_width+xspacing) ...
    yspacing*1.5+(nRows-RowIdx)*(axis_height+yspacing) axis_width axis_height],...
    holded,extravars{:});

if nargout>0
    h = temph;
end