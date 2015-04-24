function BiasValue = TestBias_multiplates(t_data, BiasCutoff, plotting, valvars)
% t_data:     need variables 'Well' or 'Row'/'Column'; all other
%                         columns will be checked against.
%  
% BiasCutoff: [minimal drop cutoff      p-value cutoff];
%                         default = [.1 .01]
%                     can be a 3x2 matrix for each condition: 
%                             edge, column, row
                            
                            
if ~exist('BiasCutoff','var') || isempty(BiasCutoff)
    BiasCutoff = [1;1;1]*[.1 .01];
elseif all(size(BiasCutoff)==[1 2])
    BiasCutoff = [1;1;1]*BiasCutoff;
elseif ~all(size(BiasCutoff)==[3 2])
    error('Wrong size for BiasCutoff input, needs to be 1x2 or 3x2')
end


if ~exist('plotting','var') || isempty(plotting)
    plotting = .5;
end

if ~exist('valvars','var') || isempty(valvars)
    valvars = 'Cellcount';
end

plates = unique(t_data.Barcode);

BiasValue = zeros(length(plates), 3);

for ip = 1:length(plates)
    t_plate = t_data(t_data.Barcode==plates(ip),intersect(varnames(t_data), ...
        [{'Barcode' 'Column' 'Row' 'Well'} valvars]));
    
    bias_res = cell(1,3);
    [biased, bias_res{:}] = TestPlateBias(t_plate, BiasCutoff, plotting);
    
    for i=find(biased)
        BiasValue(ip,i) = min(BiasValue(ip,i), min(bias_res{i}(...
            bias_res{i}(:,1)<-BiasCutoff(i,1) & bias_res{i}(:,3)<BiasCutoff(i,2),1)));
    end
end

