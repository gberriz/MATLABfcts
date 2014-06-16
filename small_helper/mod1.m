function x = mod1(x,y)
% mod1(x,y) = mod(x-1,y)+1
%   module1 return a number within [1, y+1[ instead of [0, y[ 
%       useful for indices
x = mod(x-1,y)+1;