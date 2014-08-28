function t_design = DrugDesignToTable(Design1, Perturbations, DrugOrder)
% t_design = DrugDesignToTable(Design, Perturbations, DrugOrder)
%   Convert a Design structure (not an array) with standard fields for
%   treatment into a design table. Order of Perturbations can be enforced
%   as well as Priority for the assignment of drugs in DrugName1/2/..
%
%
%   Design :    array of design structures with the following fields:
%                   - plate_dims (plate dimension)
%                   - treated_wells (wells treated with DMSO/drug/perturbation)
%                   - Drugs (structure with DrugName
%                       and layout - concentration given in uM)
%                   - Perturbations (structure with Name
%                       and layout - numeric array)
%
%   Perturbations:  list of Perturbations found in Design
%
%   DrugOrder : list of drugs to prioritize for storing in DrugName,
%                   DrugName2/... ; useful for comparison in designs with
%                   multiple drugs per well.
%
%   t_design :  table with the following columns:
%                   - DrugName (and HMSLid)
%                   - Well
%                   - Conc (for concentration in uM)
%               and annotation columns:
%                   - DMSO (added to the treated_wells)
%               and optional columns (stored as 'Perturbations')
%                   - pert_type (e.g. trt_cp, vehicle_ctl, ...)
%                   - SeedingNumber
%                   - other perturbations (e.g. EGF/...)
%                   - DrugName2/3/... and Conc2/3/... for multiple drugs per well
%


%%
[rows, cols] = find(Design1.treated_wells);
Well = ConvertRowColToWells(rows, cols);

DrugNames = {Design1.Drugs.DrugName};
if exist('DrugOrder','var')
    [temp, order] = ismember(DrugOrder, DrugNames);
    order = [order(temp) find(~ismember(DrugNames, DrugOrder))];
    DrugNames = DrugNames(order);
else
    order = 1:length(DrugNames);
end

PertNames = {Design1.Perturbations.Name};
if exist('Perturbations','var') && ~isempty(Perturbations)
    PertNames = intersect(Perturbations, PertNames, 'stable');
    if isempty(PertNames)
        warnprintf('No perturbation found in Design that matches ''Perturbations''');
    end
end

if isfield(Design1.Drugs,'HMSLid')
    t_HMSLids = table({Design1.Drugs.DrugName}', {Design1.Drugs.HMSLid}', ...
        'VariableNames', {'DrugName' 'HMSLid'});
end

%%
Ndrugs = 1;

DrugConc = reshape([Design1.Drugs(order).layout], ...
    [Design1.plate_dims length(DrugNames)]);
if any(any(sum(DrugConc>0,3)>Ndrugs))
    Ndrugs = max(max(sum(DrugConc>0,3)));
    warnprintf('some wells have %i drugs, additional columns in output', ...
        Ndrugs)
end

Conc = NaN(height(t_data),1);
DrugName = repmat({''}, height(t_data),1);
for iAD = 2:Ndrugs
    eval(sprintf('Conc%i = Conc; DrugName%i = DrugName;', ...
        iAD, iAD))
end

