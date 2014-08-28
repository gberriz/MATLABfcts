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
for iTf = 1:length(Trtfiles)
    assert(exist(fullfile(folder, Trtfiles{iTf}), 'file')>0, ...
        'Treatment file %s missing in folder %s', Trtfiles{iTf}, folder)
    [~,~,ext] = fileparts(Trtfiles{iTf});
    
    if strcmp(ext,'.mat')
        temp = load(fullfile(folder, Trtfiles{iTf}));
        Designs{iTf} = temp.Design;
    elseif strcmp(ext,'.hpdd')
        Designs{iTf} = hpdd_importer(hpdd_filename);
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

Conc = NaN(height(t_data),1);
DrugName = repmat({''}, height(t_data),1);
HMSLid = repmat({''}, height(t_data),1);

% this is not the optimal way of storing multiple drugs because of the
% hierarcy between DrugName and Conc as well as the redudancy and possible
% swapping between Drug1 and Drug2 ; it makes matching between condition hard
for iAD = 2:Ndrugs
    eval(sprintf('Conc%i = Conc; DrugName%i = DrugName; HMSLid%i = HMSLid;', ...
        iAD, iAD, iAD))
end

if exist('fields','var')
    assert(all(ismember(fields, Perturbations)), ...
        'Not all ''fields'' found as perturbations in the design files')
else
    fields = Perturbations;
end
datafields = cell(1, length(fields));


%%
% if exist('fields','var')
%     warning('need better implementation of the design file')
%     
%     temp = load(fullfile(folder, Trtfiles{1}));
%     designs = temp.design;
%     idx = t_data.Treatmentfile=='-';
%     
%     LayoutIdx = sub2ind(size(designs{iD}(1).layout), t_data.Row(idx), t_data.Column(idx));
%     
%     for iF = 1:length(fields)
%         temp = designs{1}(1).(fields{iF})(LayoutIdx);
%         datafields{iF}(idx) = temp;
%     end
% end

for iTf = 1:length(Trtfiles)        
    
    Design = Designs{iTf};
    
    
    
    
    if ~isempty(strfind(Trtfiles{iTf}, '.tsv')) || ...
            ~isempty(strfind(Trtfiles{iTf}, '.txt'))
        % case of a .tsv file
        idx = t_data.Treatmentfile==Trtfiles{iTf};
        
        t_trt = TableToCategorical(tsv2table(fullfile(folder, Trtfiles{iTf})),0);
        
        [temp, idx2] = outerjoin(t_data(idx,:), t_trt, 'keys', {'Well'}, 'rightvariables', ...
            intersect(varnames(t_trt), {'HMSLid' 'DrugName' 'Conc'}), 'Type','left');
        temp = temp(sortidx(idx2),:);
        
        DrugName(idx) = cellstr(temp.DrugName);
        Conc(idx) = temp.Conc;
        
    elseif ~isempty(strfind(Trtfiles{iTf}, '.mat'))
        % case of a MATLAB file with a structure 'design'
        
        temp = load(fullfile(folder, Trtfiles{iTf}));
        designs = temp.design;
        
        DesignNumbers = setdiff(unique(t_data.DesignNumber),0);
        DesignNumbers = DesignNumbers(~isnan(DesignNumbers));
        for iD = DesignNumbers'
            idx = t_data.Treatmentfile==Trtfiles{iTf} & t_data.DesignNumber==iD;
            
            nDrugs = length(designs{iD});
            DrugNames = [{''}; ReducName({designs{iD}.name}, ',.')'];
            DrugConc = reshape([designs{iD}.layout], ...
                [size(designs{iD}(1).layout) nDrugs]);
            DrugIdx = sum((DrugConc>0).* ...
                repmat(reshape(1:nDrugs,1,1,[]),[size(designs{iD}(1).layout), 1]),3);
            DrugConc = sum(DrugConc,3);
            LayoutIdx = sub2ind(size(designs{iD}(1).layout), t_data.Row(idx), t_data.Column(idx));
            
            DrugName(idx) = DrugNames(DrugIdx(LayoutIdx)+1);
            Conc(idx) = DrugConc(LayoutIdx);
            
            if exist('fields','var')
                for iF = 1:length(fields)
                    temp = designs{iD}(1).(fields{iF})(LayoutIdx);
                    datafields{iF}(idx) = temp;
                end
            end
            
        end
    end
end
assert(all(strcmp(DrugName(Conc==0),'') | strcmp(DrugName(Conc==0),'DMSO')))

pert_type = repmat({''}, height(t_data),1);
pert_type(Conc==0 & ismember(t_data.Treatmentfile,Trtfiles)) = {'ctl_vehicle'};
pert_type(Conc>0 & ismember(t_data.Treatmentfile,Trtfiles)) = {'trt_cp'};


t_annotated = [t_data table(DrugName, Conc, pert_type) ];
if exist('fields','var')
    datafields = cellfun2(@ToColumn, datafields);
    t_annotated = [t_annotated table(datafields{:}, 'variablenames', fields)];
end

t_annotated = TableToCategorical(t_annotated,0);



