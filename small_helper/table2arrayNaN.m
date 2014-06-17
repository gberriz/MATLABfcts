function a = table2arrayNaN(t)

IsNumeric = false(size(t,2));
for i=1:size(t,2)
    IsNumeric(i) = isnumeric(t.(t.Properties.VariableNames{i}));
end
a = NaN(size(t));
a(:,IsNumeric) = table2array(t(:,IsNumeric));
