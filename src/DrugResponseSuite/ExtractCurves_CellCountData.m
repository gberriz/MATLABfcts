function t_fits = ExtractCurves_CellCountData(t_data, keys, pcutoff)
% t_fits = ExtractCurves_CellCountData(t_data, keys, pcutoff)
%   Sigmoidal fit on the drug response data (expect concentration in uM) to
%   extract the following parameters:
%       - IC50
%       - Einf (effect -> = 0 if no drug response)
%       - Hill
%       - Area (not normalized)
%       - EC50
%   If the relative growth is measured, additional parameters are found:
%       - GI50
%       - GIinf
%   All the outputs are saved in a table with annotations including drug
%   concentrations and initial values.
%
%   keys are used for aggregation of the data; default are : CellLine,
%   DrugName, Time, SeedingNumber, Date. 
%
%   pcutoff is the cutoff for the p-value of a F-test against a flat line.
%


%%
if exist('keys','var')
    keys = intersect(t_data.Properties.VariableNames, ...
        [{'CellLine' 'DrugName' 'Time' 'Date'} keys]);    
else
    keys = intersect(t_data.Properties.VariableNames, ...
        {'CellLine' 'DrugName' 'Time' 'SeedingNumber' 'Date'});
end

if ~isempty(setdiff(intersect(varnames(t_data), ...
        strcat('DrugName', cellfun(@(x) {num2str(x)}, num2cell(2:9)))), ...
        keys))
    warnprintf('Mutiple drugs found in the data, but only the first one is a key')
end

if ~exist('pcutoff','var')
    pcutoff = .05;
end
fitopt.pcutoff = pcutoff;


t_keys = unique(t_data(:,keys));

DoGI50 = ismember('RelGrowth', t_data.Properties.VariableNames);

t_fits = table;
for ik = 1:height(t_keys)
    loop_waitbar(ik, height(t_keys))
    %%
    subt = t_data(eqtable(t_keys(ik,:), t_data(:,keys)),:);
    
    if all(subt.Conc<1e-3) || all(subt.Conc>1e3)
        warnprintf('Concentrations are expected in uM; fitopt have constraints')
    end
    
    [IC50, Hill, Einf, Area, r2, EC50, fit] = ...
        ICcurve_fit(subt.Conc, subt.RelCellCnt, 'IC50', fitopt);
    
    t_temp = [t_keys(ik,:) table(IC50, Hill, Einf, Area, r2, EC50) ...
        table({fit}, {subt.Conc'}, {subt.RelCellCnt'}, 'VariableNames', ...
        {'fit' 'Conc' 'RelCellCnt'})];
    
    if DoGI50
        [GI50, ~, GIinf, GIArea, GI_r2, ~, GI_fit] = ...
            ICcurve_fit(subt.Conc, subt.RelGrowth, 'GI50', fitopt); 
        
        t_temp = [t_temp table(GI50, GIinf, GIArea, GI_r2) ...
        table({GI_fit}, {subt.RelGrowth'}, 'VariableNames', {'GI_fit' 'RelGrowth'})];
    
    end
    
    t_fits = [t_fits; t_temp];

end