function table2tsv(a,filename)

cell2tsv(filename, [a.Properties.VariableNames;
    table2cell(a)])