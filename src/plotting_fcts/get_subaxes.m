function h = get_subaxes(nRows, nCols, RowIdx, ColIdx, holded, varargin)
% h = get_subaxes(Nrows, Ncols, RowIdx, ColIdx, holded, varargin)
%   generate new axes such that it will fit (Nrows, Ncols) axes.
%
%   (Nrows, Ncols) are the number of rows and columns
%   (RowIdx, ColIdx) are the indexes for row and column ; if ColIdx=0, will
%       behave as the MATLAB subplot (filling up each row one after
%       another)
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

if ~exist('ColIdx','var') || isempty(ColIdx) || ColIdx==0
    ColIdx = mod(RowIdx-1,nCols)+1;
    RowIdx = ceil(RowIdx/nCols);
end

p = inputParser;
p.KeepUnmatched = true;
addParameter(p, 'xspacing', .05, @isnumeric)
addParameter(p, 'yspacing', .08, @isnumeric)
addParameter(p, 'xshift', NaN, @isnumeric)
addParameter(p, 'yshift', NaN, @isnumeric)

parse(p,varargin{:})
extra = p.Unmatched;
p = p.Results;

extravars = {};
extrafields = fieldnames(extra);
for i=1:length(extrafields)
    extravars{2*i-1} = extrafields{i};
    extravars{2*i} = extra.(extrafields{i});
end

if isnan(p.xshift)
    p.xshift = 1.5*p.xspacing;
end
if isnan(p.yshift)
    p.yshift = 1.2*p.yspacing;
end

axis_width = (1-((nCols-.8)*p.xspacing+p.xshift))/nCols;
axis_height = (1-((nRows-.3)*p.yspacing+p.yshift))/nRows;

temph = get_newaxes([p.xshift+(ColIdx-1)*(axis_width+p.xspacing) ...
    p.yshift+(nRows-RowIdx)*(axis_height+p.yspacing) axis_width axis_height],...
    holded,extravars{:});

if nargout>0
    h = temph;
end