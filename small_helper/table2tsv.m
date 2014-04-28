function table2tsv(a,filename)
% table2tsv(a,filename)

cell2tsv(filename, table2cellstr(a) );