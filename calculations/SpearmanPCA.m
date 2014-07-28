function [coeff, scores, explained, latent] = SpearmanPCA(data, k)

r = corr(data,'type','spearman');

coeff = zeros(size(data,2));
latent = zeros(size(data,2),1);
explained = zeros(size(data,2),1);

idx = ~all(isnan(r));

if exist('k','var')
    assert(k<sum(idx))
    [~,l,c] = svds(r(idx,idx), k);
    latent(1:k) = diag(l);
    
    explained = 100*latent/sum(latent);
    
    % Enforce a sign convention on the coefficients -- the largest element in each
    % column will have a positive sign.
    [p,d] = size(c);
    [~,maxind] = max(abs(c),[],1);
    colsign = sign(c(maxind + (0:p:(d-1)*p)));
    coeff(idx,1:k) = bsxfun(@times,c,colsign);
    
else
    if sum(idx)>5e3
        warning('Option iwth limited number of components will be faster')
    end
    [coeff(idx,1:sum(idx)), latent(1:sum(idx)), explained(1:sum(idx))] = pcacov(single(r(idx,idx)));
end

scores = (data-repmat(mean(data),size(data,1),1))*coeff;




