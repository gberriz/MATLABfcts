
function t_Tmean = TimeMerge_CellCountData(t_mean, NTimePlates, cond_inkeys, numericfields)
% t_Tmean = TimeMerge_CellCountData(t_mean, NTimePlates, cond_inkeys, numericfields)

%% assign and control the variables
if exist('cond_inkeys','var') && ~isempty(cond_inkeys)
    cond_keys = unique([{'CellLine' 'Time' 'DrugName' 'Conc'} cond_inkeys]);
else
    cond_keys = {'CellLine' 'Time' 'DrugName' 'Conc'};
    cond_inkeys = {};
end

Relvars = intersect({'RelCellCnt' 'RelGrowth' 'nRelGrowth'}, varnames(t_mean));
labelfields = {'DesignNumber' 'Barcode' ...
    'Date' 'Row' 'Column' 'Well' 'TreatmentFile' 'Replicate'}; % remove all technical replicate info
if ~exist('numericfields','var')
    numericfields = setdiff(t_mean.Properties.VariableNames( ...
        all(cellfun(@isnumeric, table2cell(t_mean(1:min(40,end),:))))), [cond_keys labelfields Relvars]);
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

Times = unique(t_mean.Time);
Time = NaN*t_mean.Time;
% find the best plate grouping
mindiff = NaN(NTimePlates,1);
for i=1:NTimePlates
    mindiff(i) = sum(sum(diff(reshape(Times(i:(i+...
        NTimePlates*floor((length(Times)-i)/NTimePlates)-1)),NTimePlates,[]))));
end

for iT = 0:ceil(length(Times)/NTimePlates)
    t = Times(max(1,argmin(mindiff)+NTimePlates*(iT-1)):...
        min(end,argmin(mindiff)+NTimePlates*iT-1));
    Time(ismember(t_mean.Time, t)) = mean(t);
end
t_Tmean = [table(Time) t_mean(:,setdiff(varnames(t_mean),'Time','stable'))];

t_Tmean = collapse(t_Tmean, ...
    @mean, 'keyvars', setdiff(varnames(t_Tmean), [labelfields Relvars numericfields],'stable'), ...
    'valvars', [Relvars numericfields]);

t_Tmean = sortrows(t_Tmean, cond_keys);
