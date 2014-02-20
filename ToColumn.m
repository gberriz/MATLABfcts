function col = ToColumn(vect)
assert(any(size(vect)==1))
if isrow(vect);col=vect';else col=vect;end