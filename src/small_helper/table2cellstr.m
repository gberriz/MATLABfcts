function cstr = table2cellstr(t_in, headers)
% cstr = table2cellstr(t_in, headers)

t_in = TableToString(t_in, 0);

if ~exist('headers','var') || headers
    cstr = [t_in.Properties.VariableNames;    table2cell(t_in)];
else
    cstr = table2cell(t_in);
end
