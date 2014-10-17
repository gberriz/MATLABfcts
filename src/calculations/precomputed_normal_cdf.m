
function cdfprob = precomputed_normal_cdf(x,m,s)

switch nargin
    case 2
        x = x-m;
    case 3
        x = ((x-m)/s);
end

global pre_comp_normal_cdf
if isempty(pre_comp_normal_cdf)
    step = 2e-5;
    step2 = step/1e3;
    pre_comp_normal_cdf.y = unique([0 (step/2) step:step:(1-step) (1-(step/2)) 1 ...
        .001+( (-step-step2/2):(step2):(step+step2/2)) ...
        .01+( (-step-step2/2):(step2):(step+step2/2)) ...
        .05+( (-step-step2/2):(step2):(step+step2/2)) ...
        .95+( (-step-step2/2):(step2):(step+step2/2)) ...
        .99+( (-step-step2/2):(step2):(step+step2/2)) ...
        .999+( (-step-step2/2):(step2):(step+step2/2)) ], 'sorted');
    pre_comp_normal_cdf.x = icdf('normal',pre_comp_normal_cdf.y,0,1);
end
% cdfprob = pre_comp_normal_cdf.y(argmin(abs( x-pre_comp_normal_cdf.x)));
idx = find( pre_comp_normal_cdf.x>x,1,'first');
if isempty(idx), cdfprob=1; return; end
cdfprob = pre_comp_normal_cdf.y(idx-1) + diff(pre_comp_normal_cdf.y(idx-[1 0]))*...
    (x-pre_comp_normal_cdf.x(idx-1))/diff(pre_comp_normal_cdf.x(idx-[1 0]));
