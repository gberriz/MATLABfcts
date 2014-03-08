function fhandle = get_newfigure(number, position, varargin)
% fhandle = get_newfigure(number, position, varargin)
%   get a clean new figure

if ~exist('position','var')
    position = [100 80 600 450];
end

if nargout==0
    figure(number)    
else
    fhandle = figure(number);
end
clf
set(gcf, 'position', position, varargin{:});