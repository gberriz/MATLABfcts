function t_out = splitHMSLid(t_in, varname)

if~exist('varname','var')
    varname = 'DrugName';
end
    
DrugName = t_in.(varname);
iscat = iscategorical(DrugName);
if iscat    
    DrugName = cellstr(DrugName);
    DrugName(strcmp(DrugName, '<undefined>')) = {''};
end

HMSLid = cellcell2cellstr( ...
    cellfun2(@(x) regexp(x, '(HMSL\d*)', 'tokens'), DrugName));

%%
DrugName = cellcell2cellstr(cellfun2(@(x) regexp(x, 'HMSL\d*', 'split'), DrugName));
DrugName = cellfun2(@strcat,DrugName);

if iscat
    HMSLid = categorical(HMSLid);
    HMSLid(isundefined(HMSLid)) = '-';
    DrugName = categorical(DrugName);
    DrugName(isundefined(DrugName)) = '-';
end

idx = find(strcmp(t_in.Properties.VariableNames, varname));

t_out = [t_in(:,1:(idx-1)) table(DrugName, HMSLid) t_in(:,(idx+1):end)];

