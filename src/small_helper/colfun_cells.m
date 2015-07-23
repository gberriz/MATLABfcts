function b = colfun_cells(fct, a)
% b = colfun_cells(fct, a)
%   apply the function fct on each column of the cell array a to yield a
%   row cell vector b

b = cell(1,size(a,2));
for i=1:length(b)
    b{i} = fct(a(:,i));
end
