function row = ToRow(vect)
assert(any(size(vect)==1))
if iscolumn(vect);row=vect';else row=vect;end