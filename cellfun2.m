function output = cellfun2(fhandle, varargin)
output = cellfun(fhandle, varargin{:}, 'uniformoutput', false);