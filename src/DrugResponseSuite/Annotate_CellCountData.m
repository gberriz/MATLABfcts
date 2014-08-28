function t_annotated = AddTreatmentData(t_data, folder, fields)
%
%
%   annotate the data using the treatment files
%   adding the fields: 
%       - DrugName    (assuming only one drug per well)
%       - Conc
%       - any additional fields as given in input varaible 'fields' (by
%           default all 'Perturbations' in the array of structures Design)
%
%   variable folder is to specify where the TreatmentFiles are stored.
%


if ~exist('folder','var')
    folder = '';
end

%%
Trtfiles = setdiff(cellstr(unique(t_data.TreatmentFile)),'-');

Ndrugs = 1;
Perturbations = {};
% check for the treatment files
for iTf = 1:length(Trtfiles)
    assert(exist(fullfile(folder, Trtfiles{iTf}), 'file')>0, ...
        'Treatment file %s missing in folder %s', Trtfiles{iTf}, folder)
    [~,~,ext] = fileparts(Trtfiles{iTf});
    
    if strcmp(ext,'.mat')
        design_vars = load(fullfile(folder, Trtfiles{iTf}));
        Design = design_vars.Design;
    elseif strcmp(ext,'.hpdd')
        Design = hpdd_importer(hpdd_filename);
    else
        tsv2table(fullfile(folder, Trtfiles{iTf}));
        continue
    end
    
    for iD=1:size(Design)
        % check for the number of drugs in the same well
        DrugConc = reshape([Design(iD).Drugs.layout], [Design(iD).plate_dims ...
            length(Design(iD).Drugs)]);
        if any(any(sum(DrugConc>0,3)>Ndrugs))
            Ndrugs = max(max(sum(DrugConc>0,3)));
            warning('some wells have %i drugs, additional columns in output', ...
                Ndrugs)            
        end
        
        if isfield(Design(iD), 'Perturbations')
            Perturbations = unique([Perturbations {Design(iD).Perturbations.Name}], 'stable');
        end
    end
    
end
    
Conc = NaN(height(t_data),1);
DrugName = repmat({''}, height(t_data),1);

% this is not the optimal way of storing multiple drugs because of the
% hierarcy between DrugName and Conc as well as the redudancy and possible
% swapping between Drug1 and Drug2 ; it makes matching between condition hard

for iAD = 2:Ndrugs
    eval(sprintf('Conc%i = Conc; DrugName%i = DrugName;', iAD, iAD))
end

if exist('fields','var')
    datafields = cell(1, length(fields));
else
    datafields = Perturbations;
end


%%
if exist('fields','var')
    warning('need better implementation of the design file')
    
    design_vars = load(fullfile(folder, Trtfiles{1}));
    designs = design_vars.design;
    idx = t_data.Treatmentfile=='-';
    
    LayoutIdx = sub2ind(size(designs{iD}(1).layout), t_data.Row(idx), t_data.Column(idx));
    
    for iF = 1:length(fields)
        temp = designs{1}(1).(fields{iF})(LayoutIdx);
        datafields{iF}(idx) = temp;
    end
end

for iTf = 1:length(Trtfiles)
        
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
        
        design_vars = load(fullfile(folder, Trtfiles{iTf}));
        designs = design_vars.design;
        
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



