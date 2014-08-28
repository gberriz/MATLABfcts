function Design = TextDesignFile_importer(filename)
% Design = TextDesignFile_importer(filename)
%   Read a tsv file and convert the information in an array of Design structure
%   with standard fields for treatment. Assume only one plate therefore a
%   single entry per well.
%
%   filename :  path and file name to a treatment saved as tsv file with
%               the following columns:
%                   - DrugName (and HMSLid)
%                   - Well
%                   - Conc (for concentration in uM)
%               and annotation columns:
%                   - DMSO (added to the treated_wells)
%               and optional columns (stored as 'Perturbations')
%                   - TrtType (e.g. rtr_cp, vehicle_ctl, ...)
%                   - SeedingNumber
%                   - other perturbations (e.g. EGF/...)
%
%   Design :    array of design structures with the following fields:
%                   - plate_dims (plate dimension)
%                   - treated_wells (wells treated with DMSO)
%                   - well_volume (in uL)
%                   - Drugs (structure with DrugName
%                       and layout - concentration given in uM)
%                   - Perturbations (structure with Name
%                       and layout - numeric array)
%


assert(exist(filename, 'file')>0, 'Design file %s missing', filename)

t_design = tsv2table(filename);

assert(all(isvariable(t_design, {'DrugName' 'Well' 'Conc'})), ...
    'Need at least ''DrugName'' ''Well'' ''Conc'' as headers')

%% check the wells and find the size of the plate
wells = t_design.Well;
assert(length(unique(wells)) == length(wells), ...
    'Multiple values for the same well found')


t_raw.Well(well_idx) = strcat(cellcell2cellstr(regexp(t_raw.Well(well_idx),'^(\w)\d$','tokens')), ...
    '0',cellcell2cellstr(regexp(t_raw.Well(well_idx),'^\w(\d)$','tokens')));

%%

Design =  struct('plate_dims', plate_dims, 'treated_wells', treated_wells, ...
    'well_volume', well_volume, 'Drugs', Drugs);

