function b = table2cellstr(a)

tic
for i=1:size(a,2)
    if iscategorical(a.(a.Properties.VariableNames{i})(1))
        a.(a.Properties.VariableNames{i}) = cellstr(a.(a.Properties.VariableNames{i}));
    elseif isnumeric(a.(a.Properties.VariableNames{i})(1))
        a.(a.Properties.VariableNames{i}) = cellfun_brdcast(@num2str,...
            num2cell(a.(a.Properties.VariableNames{i})));
    end
end

b = [a.Properties.VariableNames;    table2cell(a)];
toc