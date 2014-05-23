function SEM_value = SEM(values)

if length(values)>2
    SEM_value = std(values)/sqrt(length(values));
else
    SEM_value = abs(diff(values)/2);
end