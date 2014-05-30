function b = table2cellstr(a)
    if height(a) > 0
        row1 = a(1, :);
        for j = 1:width(row1)
            v = row1.(j);
            if iscategorical(v) || ischar(v)
                a.(j) = cellstr(a.(j));
            elseif isnumeric(v)
                a.(j) = cellstr(num2str(a.(j)));
            end
        end
    end

    b = [a.Properties.VariableNames; table2cell(a)];
end