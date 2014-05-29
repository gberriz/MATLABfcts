function ExportProtocol_D300( base_filename, experiment )
%EXPORTPROTOCOL_D300 Write D300 protocol to '.hpdd' XML file.
%
%   EXPORTPROTOCOL_D300( filename, experiment) writes the plate maps
%   defined in experiment to filename + '.hpdd'.
%
%   experiment is a struct array (name, plate_map).
%
%   Units:
%
%      experiment(m).plate_map(n).nominal_conc  - millimolar
%      experiment(m).plate_map(n).well_volume   - microliters
%      experiment(m).plate_map(n).layout        - micromolar
%
%      (These values will undergo unit conversion as needed to meet the
%      requirements of the hpdd file format.)

document = com.mathworks.xml.XMLUtils.createDocument('Protocol');
protocol = document.getDocumentElement;

% Insert a comment to help document how the file was created.
protocol.appendChild(document.createTextNode(sprintf('\n   ')));
protocol.appendChild(document.createComment('Created by ExportProtocol_D300.m'));

% Add protocol-global settings.
create_text_children(protocol, ...
    {
    'Version'                   int2str(2);
    'VolumeUnit'                'nL';
    'ConcentrationUnit'         'µM';
    'MolarityConcentrationUnit' 'µM';
    'MassConcentrationUnit'     'ng_mL';
    'ShakePerFluid'             logical2str(0);
    'ShakePlateDuration'        int2str(5);
    'ShakePerWell'              logical2str(0);
    'ShakeThresholdVolume'      int2str(100);
    'BackfillOrder'             'LastPriority';
    'BackfillNoDispense'        logical2str(0);
    } ...
);

% Add Fluids list.
fluids = document.createElement('Fluids');
protocol.appendChild(fluids);
fluid_data = get_substances(experiment);
% Map from substance name to string id for XML output.
fluid_ids = containers.Map;
for fluid_num = 1:length(fluid_data)
    fluid = document.createElement('Fluid');
    fluids.appendChild(fluid);
    % Use 0-based numbering for fluid IDs.
    id = int2str(fluid_num - 1);
    % The utility of RelatedID is not clear to me, but it was always -1 in
    % the sample file. -JLM
    related_id = int2str(-1);
    % Add to name->id map for later.
    fluid_ids(fluid_data(fluid_num).name) = id;
    fluid.setAttribute('ID', id);
    fluid.setAttribute('RelatedID', related_id);
    % Convert concentration from millimolar to micromolar.
    conc_micromolar = fluid_data(fluid_num).nominal_conc * 1e3;
    % Create fluid properties.
    create_text_children(fluid, ...
        {
        'Name'              fluid_data(fluid_num).name;
        'Concentration'     double2str(conc_micromolar);
        % It's not clear whether this is meant as the unit to use for
        % display within the D300 software, or something else, so we'll
        % just play it safe and use the same units we used up top.
        'ConcentrationUnit' 'µM';
        } ...
    );
end

% Add Plates list.
plates = document.createElement('Plates');
protocol.appendChild(plates);
% NOTE What are the other Modes and do we need to support them?
create_text_children(plates, {'Mode' 'Concentration'});
for plate_num = 1:length(experiment)
    % Verify that all plate maps are 24x16 (384) wells. If we have to
    % support different plate types in the future this would need to be
    % changed.
    if ~all(cell2mat(cellfun2(@(x) all(size(x) == [16 24]), ...
                              {experiment(plate_num).plate_map.layout})))
        me = MException('ExportProtocol_D300:plate_size_mismatch', ...
                        'Plate %d is not 24x16=384 wells', plate_num);
        throw(me);
    end
    % TODO check that plate_map.well_volume values are identical
    plate = document.createElement('Plate');
    plates.appendChild(plate);
    % Convert volume from microliters to nanoliters.
    volume_nanoliters = experiment(plate_num).plate_map(1).well_volume * 1e3;
    create_text_children(plate, ...
        {
        'PlateType'   'Default384';
        'Rows'        int2str(16);
        'Cols'        int2str(24);
        'Name'        experiment(plate_num).name;
        'AssayVolume' double2str(volume_nanoliters);
        'DMSOLimit'   double2str(0.02);
        'DontShake'   logical2str(0);
        } ...
    );
    wells = document.createElement('Wells');
    plate.appendChild(wells);
    for row = 1:16
        for column = 1:24
            well = document.createElement('Well');
            well.setAttribute('Row', int2str(row - 1));
            well.setAttribute('Col', int2str(column - 1));
            plate_map = experiment(plate_num).plate_map;
            for map_num = 1:length(plate_map)
                % Values are already in the right units, micromolar.
                conc = plate_map(map_num).layout(row, column);
                if conc > 0
                    id = fluid_ids(plate_map(map_num).name);
                    fluid = document.createElement('Fluid');
                    fluid.setAttribute('ID', id);
                    fluid.setTextContent(double2str(conc));
                    well.appendChild(fluid);
                end
            end
            if well.hasChildNodes
                wells.appendChild(well);
            end
        end
    end
    % FIXME what is this Randomize empty element for?
    plate.appendChild(document.createElement('Randomize'));
end

% TODO Allow specification of wells to not backfill.
backfills = document.createElement('Backfills');
protocol.appendChild(backfills);
for plate_num = 1:length(experiment)
    backfill = document.createElement('Backfill');
    backfills.appendChild(backfill);
    backfill.setAttribute('Type', 'ToMaxVolume');
    wells = document.createElement('Wells');
    backfill.appendChild(wells);
    plate_id = int2str(plate_num - 1);
    for row = 1:16
        for column = 1:24
            well = document.createElement('Well');
            wells.appendChild(well);
            well.setAttribute('P', plate_id);
            well.setAttribute('R', int2str(row - 1));
            well.setAttribute('C', int2str(column - 1));
        end
    end
end

xmlwrite([base_filename '.hpdd'], document);

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


function substances = get_substances( experiment )
% Get list of unique substances from experiment. Checks that nominal_conc
% is the same for all same-named substances.
substances = struct('name', {}, 'nominal_conc', {});
for plate_num = 1:length(experiment)
    for substance_num = 1:length(experiment(plate_num).plate_map)
        substance = experiment(plate_num).plate_map(substance_num);
        name = substance.name;
        nominal_conc = substance.nominal_conc;
        match = find(strcmp(name, {substances.name}));
        if ~isempty(match)
            if substances(match).nominal_conc ~= nominal_conc
                me = MException(...
                    'ExportProtocol_D300:concentration_mismatch', ...
                    'Substance %s has different nominal_conc values', name ...
                );
                throw(me);
            end
        else
            substances(end+1) = struct( ...
                'name', name, ...
                'nominal_conc', nominal_conc ...
            ); %#ok<AGROW>
        end
    end
end
end


function T = double2str(X)
% Format a double-precision float using enough decimal digits to guarantee
% lossless round-tripping.
T = sprintf('%.17g', X);
end


function T = logical2str(X)
% Format a logical as 'True' or 'False'.
if logical(X)
    T = 'True';
else
    T = 'False';
end
end
