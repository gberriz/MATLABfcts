function t_annotated = Annotate_CellCountData(t_data, folder, fields)
% t_annotated = Annotate_CellCountData(t_data, folder, fields)
%
%   annotate the data using the treatment files
%   adding the fields: 
%       - DrugName    (assuming only one drug per well)
%       - Conc
%       - selected additional fields as given in input varaible 'fields' (by
%           default ALL 'Perturbations' in the array of structures Design)
%
%   variable folder is to specify where the TreatmentFiles are stored.
%


if ~exist('folder','var') || isempty(folder)
    folder = '';
end

%% import all treatment files and read the designs
Trtfiles = setdiff(cellstr(unique(t_data.TreatmentFile)),'-');

Designs = cell(length(Trtfiles),1);
correct_barcode = Designs;
for iTf = 1:length(Trtfiles)
    assert(exist(fullfile(folder, Trtfiles{iTf}), 'file')>0, ...
        'Treatment file %s missing in folder %s', Trtfiles{iTf}, folder)
    [~,~,ext] = fileparts(Trtfiles{iTf});
    fprintf('\tLoading %s\n', Trtfiles{iTf})
    
    if strcmp(ext,'.mat')
        temp = load(fullfile(folder, Trtfiles{iTf}));
        Designs{iTf} = temp.Design;
    elseif strcmp(ext,'.hpdd')        
        [Designs{iTf}, correct_barcode{iTf}] = hpdd_importer(fullfile(folder, Trtfiles{iTf}));
        % because of redundant plates, barcodes have to be reassigned        
    else
        try % assume a tsv file
            Designs{iTf} = TextDesignFile_importer(fullfile(folder, Trtfiles{iTf}));
        catch err
            error(['Expecting a .mat, .hpdd or a formated tab-separated ' ...
                'file as TreatmentFile\nError : %s'], err)
        end
    end
    
end

%% look up all the drugs and perturbations in the designs

Ndrugs = 1;
Perturbations = {};
DrugNames = {};
t_HMSLids = table;
allDesigns = [Designs{:}];
for iD=1:size(allDesigns)
    % check for multiple drugs in the same well
    DrugConc = reshape([allDesigns(iD).Drugs.layout], [allDesigns(iD).plate_dims ...
        length(allDesigns(iD).Drugs)]);
    DrugNames = unique([DrugNames {allDesigns(iD).Drugs.DrugName}],'stable');
    if isfield(allDesigns(iD).Drugs,'HMSLid')
        t_HMSLids = [ t_HMSLids; table({allDesigns(iD).Drugs.DrugName}', ...
            {allDesigns(iD).Drugs.HMSLid}','VariableNames', {'DrugName' 'HMSLid'})];
    end
    if any(any(sum(DrugConc>0,3)>Ndrugs))
        Ndrugs = max(max(sum(DrugConc>0,3)));
        warnprintf('some wells have %i drugs, additional columns in output', ...
            Ndrugs)
    end
    
    % store all possible perturbations
    if isfield(allDesigns(iD), 'Perturbations')
        Perturbations = unique([Perturbations {allDesigns(iD).Perturbations.Name}], 'stable');
    end
end
t_HMSLids = unique(t_HMSLids);

   
%% declare the variables


% this is not the optimal way of storing multiple drugs because of the
% hierarcy between DrugName and Conc as well as the redudancy and possible
% swapping between Drug1 and Drug2 ; it makes matching between condition hard


if exist('fields','var')
    assert(all(ismember(fields, Perturbations)), ...
        'Not all ''fields'' found as perturbations in the design files')
else
    fields = Perturbations;
end
datafields = cell(1, length(fields));


%%
t_annotated = t_data(t_data.TreatmentFile=='-',:);

%%%%% issue with the untreated plates ... needs to fill up the columns to
%%%%% merge !!



for iTf = 1:length(Trtfiles)        
    
    DesignNumbers = unique(t_data.DesignNumber(t_data.TreatmentFile==Trtfiles{iTf}));
    fprintf('Design %s:\n', Trtfiles{iTf} )
    for iDN = 1:length(DesignNumbers)        
        if ~isempty(correct_barcode{iTf})
            DNidx = correct_barcode{iTf}.DesignNumber(DesignNumbers(iDN));
            fprintf('\thpdd exp %i -> design %i\n', DesignNumbers(iDN), DNidx);
        else
            DNidx = iDN;
        end
        t_design = DrugDesignToTable(Designs{iTf}(DNidx), fields, DrugNames);        
        idx = t_data.TreatmentFile==Trtfiles{iTf} & t_data.DesignNumber==DesignNumbers(iDN);
        temp = join(t_data(idx,:), t_design, 'keys', 'Well');
        assert(height(temp)==sum(idx))
        t_annotated = [t_annotated; temp];            
    end
    
end
    
assert(height(t_annotated)==height(t_data), 'table went from %i to %i rows; check labels', ...
    height(t_data), height(t_annotated))
