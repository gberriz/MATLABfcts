% IsUnique = HasUniqueElement(x)
%   true is x contains the same value/str
%   x should be compatible with MATLAB function 'unique'
%

function IsUnique = HasUniqueElement(x)

IsUnique = length(unique(x))==1;
