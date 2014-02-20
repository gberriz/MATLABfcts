function fhandle = get_newfigure(number, holded)
% fhandle = get_newfigure(number, holded)

if nargout==0
    figure(number);
else
    fhandle = figure(number);
end
clf
if exist('holded','var') && holded
    hold on
end