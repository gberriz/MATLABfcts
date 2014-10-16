
function cdfprob = precomputed_normal_cdf(x,m,s)
if ~exist('m','var')
    m = 0;
end
if ~exist('s','var')
    s = 0;
end

global pre_comp_normal_cdf
if isempty(pre_comp_normal_cdf)
    step = 2e-5;
    pre_comp_normal_cdf.y = [0 (step/2) step:step:(1-step) (1-(step/2)) 1];
    pre_comp_normal_cdf.x = icdf('normal',pre_comp_normal_cdf.y,0,1);
end
% cdfprob = pre_comp_normal_cdf.y(argmin(abs( ((x-m)/s) -pre_comp_normal_cdf.x)));
idx = find( pre_comp_normal_cdf.x>((x-m)/s),1,'first');
if isempty(idx), cdfprob=1; return; end
cdfprob = pre_comp_normal_cdf.y(idx-1) + diff(pre_comp_normal_cdf.y(idx-[1 0]))*...
    (pre_comp_normal_cdf.x(idx)-((x-m)/s))/diff(pre_comp_normal_cdf.x(idx-[1 0]));
