function table2tsv(a,filename)
% table2tsv(a,filename)

c = table2cell(a);
c(cellfun(@iscategorical,c)) = cellfun2(@char,c(cellfun(@iscategorical,c)));

cell2tsv(filename, [a.Properties.VariableNames;
    c])