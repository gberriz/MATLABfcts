function mat_out = NormUnit(mat_in, dim)

if ~exist('dim','var')
    dim = 1;
end
dims = ones(1,length(size(mat_in)));
dims(dim) = size(mat_in,dim);
 
mat_out = mat_in ./ repmat(sum(mat_in,dim), dims);