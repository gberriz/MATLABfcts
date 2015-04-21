function out = tail(x, varargin)
% out = tail(x, varargin)
%   relies on the function 'head' from datarail

out = head(flipud(x), varargin{:});
out = flipud(out);
