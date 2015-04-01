function varargout = xyerrorbars(x, xerr, y, yerr, varargin)
% h = xyerrorbars(x, xerr, y, yerr, varargin)
%   xerr and yerr can be empty or a unique value
%   varargin: same arguments as plots
%
%   h(1) = line (x,y) handle
%   h(2) = errorbar handle
%

assert(length(x)==length(xerr) || length(xerr)<=1)
assert(length(x)==length(y))
assert(length(x)==length(yerr)  || length(yerr)<=1)


ish = ishold;

h = plot(x,y,varargin{:});
% set(h,'marker','none')
hold on

errbars = NaN(0,2);
for i=1:length(x)
    if ~isempty(xerr)
        errbars = [errbars
            NaN NaN
            x(i)+xerr(min(i,length(xerr)))*[-1;1] y(i)*[1;1]
            ];
    end
    if ~isempty(yerr)
        errbars = [errbars
            NaN NaN
            x(i)*[1;1] y(i)+yerr(min(i,length(yerr)))*[-1;1]
            ];
    end
end

h(2) = plot(errbars(:,1), errbars(:,2),varargin{:});
set(h(2),'marker','none','linewidth',.5)

if ~ish, hold off, end

if nargout>0
    varargout = {h};
end
