function c = num2cellstr(x)

c = cellfun(@num2str, num2cell(x), 'uniformoutput',0);