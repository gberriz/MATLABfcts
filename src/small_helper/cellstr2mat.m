function mat = cellstr2mat(data)

numcell = cellfun2(@str2double, data);
numcell(cellfun(@isempty,numcell)) = {NaN};
mat = cell2mat(numcell);
