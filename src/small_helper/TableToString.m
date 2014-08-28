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
    if iscategorical(t_out.(ivar{:}))
        t_out.(ivar{:}) = cellstr(t_out.(ivar{:}));
    elseif isnumeric(t_out.(ivar{:}))
        t_out.(ivar{:}) = cellfun2(@num2str, num2cell(t_out.(ivar{:})));
    elseif iscellstr(t_out.(ivar{:}))
        t_out.(ivar{:}) = t_out.(ivar{:});
    else
        temp = cell(height(t_out),1);
        for i=1:length(temp)
            temp{i} = disp2str(t_out.(ivar{:}){i});
        end
        t_out.(ivar{:}) = temp;
    end
end
