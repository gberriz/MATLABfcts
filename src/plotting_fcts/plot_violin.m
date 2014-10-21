function h = plot_violin(x, y, ybins, width, varargin)

if length(x)>1
    assert(any (length(x) == size(y)), 'x and y must have a common dimension')
    if length(x) ~= size(y,1)
        idx = find(length(x) == size(y),1,'first');
        y = permute(y,[ idx setdiff(1:length(size(y)),idx)]);
    end
    
    if ~exist('ybins','var') || isempty(ybins)
        ydiff = (max(y(:))-min(y(:)))/20;
        ybins = (min(y(:))-ydiff*3):( ydiff ):(max(y(:))+ydiff*3);
    end
    
    if ~exist('width','var') || isempty(width)
        width = .7;
    end
    
    ish = ishold;
    hold on
    for i=1:length(x)
        h = plot_violin(x(i), y(i,:), ybins, width, varargin{:});
    end
    if ~ish
        hold off
    end
    return
end

if ~exist('width','var') || isempty(width)
    width = .7;
end

if ~exist('ybins','var') || isempty(ybins)
    ydiff = (max(y)-min(y))/20;
    ybins = (min(y)-ydiff*3):( ydiff ):(max(y)+ydiff*3);
else
    ydiff = diff(ybins(1:2));
end

ydist = ksdensity(y, ybins, 'width', ydiff*1.3);
ydist = ToRow(ydist)*width/2/max(ydist);

h = plot([x+ydist NaN x-ydist], [ybins NaN ybins], varargin{:});
h(2) = plot(x, median(y), '.k');