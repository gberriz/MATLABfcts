function hpdd_exporter(hpdd_filename, Designs, t_plateinfo)
%HPDD_IMPORTER(hpdd_filename, Designs, t_plateinfo)
%   Write an hpdd file based on array of Designs and a plate barcode table.
%
%   base_pathname : path and file name to a D300 protocol file ('.hpdd'
%                       automatically appended -- do not include it here)
%   Designs :       array of design structures with the following fields:
%                       - plate_dims (plate dimension)
%                       - treated_wells (wells treated with DMSO)
%                       - well_volume (in uL)
%                       - Drugs (structure with DrugName, HMSLid, stock_conc
%                           and layout - concentration given in uM)
%   t_plateinfo :   table (or file name to a tsv table) with columns:
%                       - Barcode
%                       - TreatmentFile (which should match the name where 
%                           'Deisgn is saved')
%                       - DesignNumber

document = com.mathworks.xml.XMLUtils.createDocument('Protocol');
protocol = document.getDocumentElement;

% Insert a comment to help document how the file was created.
protocol.appendChild(document.createTextNode(sprintf('\n   ')));
protocol.appendChild(document.createComment('Created by hpdd_exporter.m'));

% Constant for the non-ascii character mu.
MICRO = char(181);

% Add protocol-global settings.
create_text_children(protocol, ...
    {
    'Version'                   int2str(2);
    'VolumeUnit'                'nL';
    'ConcentrationUnit'         [MICRO 'M'];
    'MolarityConcentrationUnit' [MICRO 'M'];
    'MassConcentrationUnit'     'ng_mL';
    'ShakePerFluid'             logical2str(true);
    'ShakePlateDuration'        int2str(5);
    'ShakePerWell'              logical2str(false);
    'ShakeThresholdVolume'      int2str(100);
    'BackfillOrder'             'LastPriority';
    'BackfillNoDispense'        logical2str(0);
    } ...
);

% Add Fluids list.
fluids = document.createElement('Fluids');
protocol.appendChild(fluids);
fluid_data = get_DesignDrugs(Designs);
% Map from drug name to string id for XML output.
fluid_ids = containers.Map;
for fluid_num = 1:length(fluid_data)
    fluid = document.createElement('Fluid');
    fluids.appendChild(fluid);
    % Use 0-based numbering for fluid IDs.
    id = int2str(fluid_num - 1);
    % The utility of RelatedID is not clear to me, but it was always -1 in
    % one sample file, and not present in another. -JLM
    related_id = int2str(-1);
    % Add to name->id map for later.
    fluid_ids(fluid_data(fluid_num).name) = id;
    fluid.setAttribute('ID', id);
    fluid.setAttribute('RelatedID', related_id);
    % Concentration shoudl be micromolar.
    conc_micromolar = fluid_data(fluid_num).stock_conc;
    % Create fluid properties.
    create_text_children(fluid, ...
        {
        'Name'              fluid_data(fluid_num).name;
        'Concentration'     double2str(conc_micromolar);
        % It's not clear whether this is meant as the unit to use for
        % display within the D300 software, or something else, so we'll
        % just play it safe and use the same units we used up top.
        'ConcentrationUnit' [MICRO 'M'];
        } ...
    );
end

% Perform some checks on the Design data structures.
for design_num = 1:length(Designs)
    % Verify that all drug layout have the standard format
    if size(unique(cell2mat(cellfun(@size,{Designs(design_num).Drugs.layout}','uniformoutput',0)), ...
            'rows'),1)>1
        me = MException('ExportProtocol_D300:drug_layout_mismatch', ...
                        'Design %d have different layout sizes for drugs', design_num);
        throw(me);
    end
    % Verify that all plate maps are standard format. If we have to
    % support different plate types in the future this would need to be
    % changed.
    if (log2(size(Designs(design_num).Drugs(1).layout,1))~=...
            round(log2(size(Designs(design_num).Drugs(1).layout,1)))) || ...
            ((size(Designs(design_num).Drugs(1).layout,2)/...
            size(Designs(design_num).Drugs(1).layout,1))~=1.5)
        me = MException('ExportProtocol_D300:plate_size_mismatch', ...
                        'Design %d is not a standard plate size', design_num);
        throw(me);
    end
    % Verify that backfill maps are 16x24 (384) wells if it exists. If we have to
    % support different plate types in the future this would need to be
    % changed.
    if isfield(Designs(design_num), 'treated_wells') && ...
            ~all(size(Designs(design_num).treated_wells) == size(Designs(design_num).Drugs(1).layout))
        me = MException('ExportProtocol_D300:backfill_size_mismatch', ...
                        'treated_wells for Design %d is not matching drug layout', design_num);
        throw(me);
    end
end

% load the table is a file name was passed
if ischar(t_plateinfo)
    t_plateinfo = tsv2table(t_plateinfo);
else
    assert(istable(t_plateinfo), ['Barcodes should be a table or a tsv file' ...
        ' with columns Barcode, DesignNumber'])
end
    
% Create Plates container.
plates = document.createElement('Plates');
protocol.appendChild(plates);
% NOTE What are the other Modes and do we need to support them?
create_text_children(plates, {'Mode' 'Concentration'});

% Create Backfills container.
backfills = document.createElement('Backfills');
protocol.appendChild(backfills);
backfill = document.createElement('Backfill');
backfills.appendChild(backfill);
backfill.setAttribute('Type', 'ToMaxVolume');
backfill_wells = document.createElement('Wells');
backfill.appendChild(backfill_wells);
% FIXME Saw a ClassID="0" attribute in here. What is this?
% FIXME For a six-plate experiment, saw two backfill elements each with half of
%   the wells in it. Why?

assert(length(setdiff(unique(t_plateinfo.TreatmentFile),'-'))==1, ...
    ['Only one treatment file can be specified in the plate info file; ' ...
    'it should correspond to .mat file where the variable ''Designs'' is saved'])

t_trt_plates = t_plateinfo(~strcmp(t_plateinfo.TreatmentFile,'-'),:);

for plate_num = 1:height(t_trt_plates)
    
    plate = document.createElement('Plate');
    plates.appendChild(plate);
    % Use barcode table's DesignNumber column as index into Design array.
    cur_design = Designs(t_trt_plates.DesignNumber(plate_num));
    plate_name = t_trt_plates.Barcode(plate_num);
    if isvariable(t_trt_plates, 'PlateShaking')
        PlateShaking = t_trt_plates.PlateShaking(plate_num);
    else
        PlateShaking = false;
    end
    % Convert volume from microliters to nanoliters.
    volume_nanoliters = cur_design.well_volume * 1e3;
    create_text_children(plate, ...
        {
        'PlateType'   sprintf('Default%i', numel(cur_design.Drugs(1).layout));
        'Rows'        int2str(size(cur_design.Drugs(1).layout,1));
        'Cols'        int2str(size(cur_design.Drugs(1).layout,2));
        'Name'        plate_name;
        'AssayVolume' double2str(volume_nanoliters);
        'DMSOLimit'   double2str(0.02);
        'DontShake'   logical2str(~PlateShaking);
    % FIXME: I came across an AqueousLimit element in a different
    % officially-generated file. What is that?
        } ...
    );
    wells = document.createElement('Wells');
    plate.appendChild(wells);
    for row = 1:size(cur_design.Drugs(1).layout,1)
        for column = 1:size(cur_design.Drugs(1).layout,2)
            % Create main well/fluid elements for drug treatments.
            well = document.createElement('Well');
            well.setAttribute('Row', int2str(row - 1));
            well.setAttribute('Col', int2str(column - 1));
            drugs = cur_design.Drugs;
            for drug_num = 1:length(drugs)
                % Values are already in the right units, micromolar.
                conc = drugs(drug_num).layout(row, column);
                if conc > 0
                    id = fluid_ids(drug2displayname(drugs(drug_num)));
                    fluid = document.createElement('Fluid');
                    fluid.setAttribute('ID', id);
                    fluid.setTextContent(double2str(conc));
                    well.appendChild(fluid);
                end
            end
            % Only add well element to document if not empty.
            if well.hasChildNodes
                wells.appendChild(well);
            end
            % Create backfill well element. If treated_wells is not present we
            % apply backfill to all wells. If it is present, we look up the
            % current well address in it to determine whether to apply backfill.
            if ~isfield(cur_design, 'treated_wells') || ...
                    cur_design.treated_wells(row,column)
                backfill_well = document.createElement('Well');
                backfill_wells.appendChild(backfill_well);
                backfill_well.setAttribute('P', int2str(plate_num - 1));
                backfill_well.setAttribute('R', int2str(row - 1));
                backfill_well.setAttribute('C', int2str(column - 1));
            end
        end
    end
    % Even without randomization, the sample file had an empty Randomize tag.
    plate.appendChild(document.createElement('Randomize'));
end

xmlwrite([hpdd_filename '.hpdd'], document);

Write_DesignTreatment_summary([hpdd_filename '_summary.tsv'], Designs, t_trt_plates)

end


function create_text_children( element, data )
% Create a list of child nodes under element, using name/text pairs from
% the cell array in data.
document = element.getOwnerDocument;
for i = 1:size(data, 1)
    child = document.createElement(data{i,1});
    element.appendChild(child);
    child.setTextContent(data{i,2});
end
end




function name = drug2displayname(drug)
% Return the 'DrugName' and 'HMSLid' fields of drug, joined with a space.
name = drug.DrugName;
if length(drug.HMSLid) > 0
    name = [name ' ' drug.HMSLid];
end
end


function T = double2str(X)
% Format a double-precision float using 9 digits of precision, the same as
% the official D300 software appears to use.
T = sprintf('%.9g', X);
end


function T = logical2str(X)
% Format a logical as 'True' or 'False'.
if logical(X)
    T = 'True';
else
    T = 'False';
end
end

