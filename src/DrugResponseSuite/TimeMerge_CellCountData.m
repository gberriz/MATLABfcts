
function t_Tmean = TimeMerge_CellCountData(t_processed, NTimePlates, plate_inkeys, cond_inkeys, numericfields)
% t_Tmean = TimeMerge_CellCountData(t_processed, NTimePlates, plate_inkeys, cond_inkeys, numericfields)

%% assign and control the variables
    plate_keys = {'CellLine'};
if exist('cond_inkeys','var') && ~isempty(plate_inkeys)
    plate_keys = unique([plate_keys plate_inkeys]);
end

cond_keys = {'DrugName' 'Conc' 'Time'};
if exist('cond_inkeys','var') && ~isempty(cond_inkeys)
    cond_keys = unique([cond_keys cond_inkeys]);
end

Relvars = intersect({'RelCellCnt' 'RelGrowth' 'nRelGrowth'}, varnames(t_processed));
labelfields = {'DesignNumber' 'Barcode' ...
    'Date' 'Row' 'Column' 'Well' 'TreatmentFile' 'Replicate'}; % remove all technical replicate info
if ~exist('numericfields','var')
    numericfields = setdiff(t_processed.Properties.VariableNames( ...
        all(cellfun(@isnumeric, table2cell(t_processed(1:min(40,end),:))))), ...
        [plate_keys cond_keys labelfields Relvars]);
    if ~isempty(numericfields)
        fprintf('\tThese numeric fields will be averaged (set as cond_inkeys to use them as key):\n');
        for i=1:length(numericfields)
            fprintf('\t - %s\n', numericfields{i});
        end
    end
end

%%

% find the number of different of plates to merge and group them based on
% the plate_keys (with Time==0)

t_plates = unique(t_processed(:,plate_keys));
t_Tmean = table;

for iP = 1:height(t_plates)
    
    subt = t_processed(eqtable(t_processed, t_plates(iP,:)),:);
    
    Times = unique(subt.Time);
    Time = NaN*subt.Time;
    % find the best plate grouping
    mindiff = NaN(NTimePlates,1);
    for i=1:NTimePlates
        mindiff(i) = sum(sum(diff(reshape(Times(i:(i+...
            NTimePlates*floor((length(Times)-i)/NTimePlates)-1)),NTimePlates,[]))));
    end
    
    for iT = 0:ceil(length(Times)/NTimePlates)
        t = Times(max(1,argmin(mindiff)+NTimePlates*(iT-1)):...
            min(end,argmin(mindiff)+NTimePlates*iT-1));
        Time(ismember(subt.Time, t)) = mean(t);
    end
    Time = round(Time,1);
    temp = [table(Time) subt(:,setdiff(varnames(subt),'Time','stable'))];
    
    t_Tmean = [t_Tmean; collapse(temp, ...
        @mean, 'keyvars', setdiff([cond_keys varnames(temp)], [labelfields Relvars numericfields],'stable'), ...
        'valvars', [Relvars numericfields])];
end

t_Tmean = sortrows(t_Tmean, plate_keys);
