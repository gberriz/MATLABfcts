function theAxis = subplot2(nRows, nCols, RowId, ColId, varargin)

plotId = ColId + nCols*(RowId-1);
if nargout>0
    theAxis = subplot(nRows, nCols, plotId, varargin{:});
else
    subplot(nRows, nCols, plotId, varargin{:});
end