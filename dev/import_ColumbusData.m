% [t_CL, Drugs, CellLines] = import_ColombusData(FilesLocation, savefile)
%
%   import the results of the D300 treated plates scanned with Columbus
%
%   all path to files and parameters have to be passed in FilesLocation as:
%             FilesLocation.folder = './results_20140307/';
%                       folder with the designs and the results (in onefile)
%             FilesLocation.data_file = '20140310_TransCenter combo.txt';
%                       result file (output from Columbus)
%             FilesLocation.DesignFiles = 'Design_20140307_TransCenter_';
%                       prefix for the treatment desing (e.g.
%                           Design_20140307_TransCenter_MDA231-2.mat)                       
%             FilesLocation.CellLines = {'HME' 'MDA231'};
%                       Cell line names (prefix for the Design file)
%             FilesLocation.ColumbusTag = {'HME' '231'};
%                       Tag used in the Columbus result file for cell lines
% 
%             FilesLocation.Replicate_tags = {'drug1' 'drug2' 'drug3'};
%                       Tag used in the Columbus result file for treatment
%                       replicates
%             FilesLocation.Untrt_tag = {'day3'};   [optional]
%                       Tag used in the Columbus result file for
%                       untreated plates (controls)
%             FilesLocation.Day0_tag = {'day0'};    [optional]
%                       Tag used in the Columbus result file for 
%                       plates fixed at the day of treatment.
%             FilesLocation.OutputFile = 'Results_20140310_TransCenter';
%             [needed in using 'savefile']
%                       File name for saving the results
%   
%  savefile = save results in a file (need


function [t_CL, Drugs, CellLines] = import_ColumbusData(FilesLocation, savefile)

if ~exist('savefile','var'),
    savefile = 1;
end
    
fields = {'folder' 'data_file' 'DesignFiles'  'CellLines' ...
    'ColumbusTag' 'Replicate_tags' 'Untrt_tag' 'Day0_tag' 'OutputFile'};
Mandatory = [1   1   1   1         1   1  0   0  savefile];

for i=1:length(fields)
    if isfield(FilesLocation, fields{i})
        eval([fields{i} '=FilesLocation.' fields{i} ';'])
    elseif Mandatory(i)
        error('Field %s needed for the import', fields{i})
    else
        eval([fields{i} '={};'])
    end
end


%%

t_raw = readtable([folder data_file],'delimiter','\t');

assert(length(unique(t_raw.NumberOfAnalyzedFields))==1)

allExp = [Replicate_tags Untrt_tag Day0_tag];


assert(length(unique(t_raw.Result))==...
    (length(CellLines)*length(allExp)));

url = regexp(t_raw.URL,'/','split');
PlateID = vertcat(url{:}); 
PlateID = cellfun(@str2num,PlateID(:,end-1));

CellLine = cell(size(t_raw.Result));
Replicate = zeros(size(t_raw.Result));
Untrt = zeros(size(t_raw.Result));
Day0 = zeros(size(t_raw.Result));

cnt = 0;
for iCL = 1:length(CellLines)
    for iR = 1:length(allExp)
        
        idx = find(strfindcell(t_raw.Result, ColumbusTag{iCL})>0 & ...
            strfindcell(t_raw.Result, allExp{iR})>0);
        assert(length(unique(PlateID(idx)))==1)
        assert(isempty(intersect(unique(PlateID(idx)), unique(PlateID(~idx)))))
        
        CellLine(idx) = CellLines(iCL);
        if ismember(allExp{iR}, Replicate_tags)
            Replicate(idx) = find(strcmp(allExp{iR}, Replicate_tags));
        end
        
        if ismember(allExp{iR}, Untrt_tag)
            Untrt(idx) = find(strcmp(allExp{iR}, Untrt_tag));
        end
        
        if ismember(allExp{iR}, Day0_tag)
            Day0(idx) = find(strcmp(allExp{iR}, Day0_tag));
        end
        
        assert(all(Replicate(idx)>0 | Untrt(idx)>0 | Day0(idx)>0))
        cnt = cnt+1;
    end
end
assert(cnt==length(unique(PlateID)))

CellLine = categorical(CellLine);
t_data = [table(CellLine,Replicate,Untrt,Day0) t_raw(:,{'Well' 'Row' 'Column' 'Nuclei_NumberOfObjects'})];
t_data.Properties.VariableNames{'Nuclei_NumberOfObjects'} = 'Cellcount';

%%
t_CL = cell(1,length(CellLines));
Drugs = t_CL;

for iCL=1:length(CellLines)
    
    t_CL{iCL} = t_data(t_data.CellLine==CellLines{iCL},:);
    load([folder DesignFiles CellLines{iCL} '-1.mat'])        
    Drugs{iCL} = {drugs_struct.name};
    
    DrugDoses = NaN(height(t_CL{iCL}),length(Drugs{iCL}));
    t_CL{iCL} = sortrows(t_CL{iCL},{'Replicate' 'Column' 'Row'});
    for iR = setdiff(unique(t_data.Replicate), 0)'
        load([folder DesignFiles CellLines{iCL} '-' num2str(iR) '.mat'])
        assert(all(strcmp(Drugs{iCL}, {drugs_struct.name})))        
        for iD=1:length(Drugs{iCL})
            DrugDoses(t_CL{iCL}.Replicate==iR,iD) = drugs_struct(iD).layout(:);
        end
    end
    
    if isfield(drugs_struct,'IsPrimary')
        Drugs{iCL} = struct('DrugName',ReducName(Drugs{iCL}), ...
            'IsPrimary',{drugs_struct.Primary},'ComboDoses',{drugs_struct.Doses},...
            'SingleDoses',{drugs_struct.SingleDoses});
    else
        Drugs{iCL} = struct('DrugName',ReducName(Drugs{iCL}), ...
            'ComboDoses',{drugs_struct.Doses},...
            'SingleDoses',{drugs_struct.SingleDoses});
    end
        
    CtrlIdx = all(DrugDoses==0,2);
    for iR = setdiff(unique(t_CL{iCL}.Replicate), 0)'
        % remove corner as control to avoid bias (only the case if more
        % than 20 controls overall (enough control remain)
        CornerIdx = find(t_CL{iCL}.Replicate==iR & ...
            ismember(t_CL{iCL}.Well,{'A1' 'A24' 'P1' 'P24'}));
        assert(length(CornerIdx)==4)
        if all(CtrlIdx(CornerIdx)==1)
            CtrlIdx(CornerIdx)=0;
        end
    end
    
    t_CL{iCL} = [t_CL{iCL} array2table([DrugDoses CtrlIdx], ...
        'VariableNames',{Drugs{iCL}.DrugName 'Ctrl'})];
        
    t_CL{iCL} = sortrows(t_CL{iCL},{'Replicate' 'Column' 'Untrt' ...
        'Day0' 'Ctrl'});
end

if savefile
    save([OutputFile '.mat'], 't_CL', 'Drugs', 'CellLines')
end
