
function t_Tmean = TimeMerge_CellCountData(t_mean, NTimePlates, cond_inkeys, numericfields)
% t_Tmean = TimeMerge_CellCountData(t_mean, NTimePlates, cond_inkeys, numericfields)

%% assign and control the variables
if exist('cond_inkeys','var') && ~isempty(cond_inkeys)
    cond_keys = unique([{'CellLine' 'Time' 'DrugName' 'Conc'} cond_inkeys]);
else
    cond_keys = {'CellLine' 'Time' 'DrugName' 'Conc'};
    cond_inkeys = {};
end

labelfields = {'pert_type' 'RelCellCnt' 'RelGrowth' 'nRelGrowth' 'DesignNumber' 'Barcode' ...
    'Untrt' 'Date' 'Row' 'Column' 'Well' 'TreatmentFile' 'Replicate'};
if ~exist('numericfields','var')
    numericfields = setdiff(t_mean.Properties.VariableNames( ...
        all(cellfun(@isnumeric, table2cell(t_mean)))), [cond_keys labelfields]);
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
Relvars = intersect({'RelCellCnt' 'RelGrowth' 'nRelGrowth'}, varnames(t_mean));
t_cond = unique(t_mean(:,setdiff(varnames(t_mean),['Time' numericfields Relvars], 'stable')));

t_Tmean = table;

% loop through the different plates
for iC = 1:height(t_cond)
    %
    loop_waitbar(iC, height(t_cond))
    t_conditions = sortrows(t_mean(eqtable(t_mean, t_cond(iC,:)),:), 'Time');
    
    % find the best plate grouping
    Times = t_conditions.Time;
    mindiff = NaN(NTimePlates,1);
    for i=1:NTimePlates
        mindiff(i) = sum(sum(diff(reshape(Times(i:(end-NTimePlates+i-1)),NTimePlates,[]))));
    end
    
    %%%%%%%%% maybe move that part out of the loop if there are
    %%%%%%%%% corresponding replicates from different plates
    
    %% collapse the replicates
    % can be done much smarter and faster by rewriting it as an array and
    % averaging it.
    for iT = 0:ceil(length(Times)/NTimePlates)
        temp = collapse(t_conditions(max(1,argmin(mindiff)+NTimePlates*(iT-1)):...
            min(end,argmin(mindiff)+NTimePlates*iT-1),:), ...
            @mean, 'keyvars', setdiff(varnames(t_mean), [Relvars 'Time' numericfields],'stable'), ...
            'valvars', [Relvars 'Time' numericfields]);
        
        t_Tmean = [t_Tmean; temp];
    end
end


t_Tmean = sortrows(t_Tmean, cond_keys);

end
