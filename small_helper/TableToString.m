function t_out = TableToString(t_in, varidx)
% TableToCategorical(t_in, varidx)

if isnumeric(varidx)
    if varidx==0
        varidx = t_in.Properties.VariableNames;
    else
        varidx = t_in.Properties.VariableNames(varidx);
    end
end

t_out = t_in;

for ivar = ToRow(varidx)
    t_out.(ivar{:}) = cellstr(t_out.(ivar{:}));
end
