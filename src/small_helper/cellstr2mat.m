function mat = strcell2mat(data)

numcell = cellfun2(@str2num, data);
numcell(cellfun(@isempty,numcell)) = {NaN};
mat = cell2mat(numcell);