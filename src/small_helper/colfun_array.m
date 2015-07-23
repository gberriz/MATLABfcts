function b = colfun_array(fct, a)
% b = colfun_array(fct, a)
%   apply the function fct on each column of the matrix a to yield a row
%   vector b

b = NaN(1,size(a,2));
for i=1:length(b)
    b(i) = fct(a(:,i));
end
