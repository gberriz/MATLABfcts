function [coeff, scores, explained, latent] = SpearmanPCA(data)

r = corr(data,'type','spearman');

coeff = zeros(size(r));
latent = zeros(size(r,1),1);
explained = zeros(size(r,1),1);

idx = ~all(isnan(r));

[coeff(idx,1:length(idx)), latent(1:length(idx)), explained(1:length(idx))] = pcacov(r(idx,idx));

scores = zeros(size(data));
scores(:,idx) = data(:,idx)*coeff(:,idx);