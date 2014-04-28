function cell2csv(filename,data,delimiter)
%  cell2csv(filename,data,delimiter)
%       default delimiter = ','


if ~exist('delimiter','var')
    delimiter = ',';
end
if delimiter==' '
    error('not working, needs to be corrected')
end

numidx = cellfun(@isnumeric,data(:));
data(numidx) = cellfun(@num2str,data(numidx),'uniformoutput',0);

data(:,1:end-1) = strcat(data(:,1:end-1),delimiter);

file = fopen(filename,'w');
for i=1:size(data,1)
    fprintf(file,[data{i,:}]);
    if i<size(data,1)
        fprintf(file,'\n');
    end
end

pause(.1)
fclose(file);