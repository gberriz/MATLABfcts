function out = cellstrcat(x, delimiter)

if ~exist('delimiter','var')
    delimiter = ' ';
end

x = x(~cellfun(@isempty, x));

if isempty(x)
    out = '';
else
    out = x{1};
    for i=1:numel(x)
        out = [out delimiter x{i}];
    end
end
