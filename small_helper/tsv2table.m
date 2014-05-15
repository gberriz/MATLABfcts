function t = tsv2table(filename)

fprintf('Importing %s\n', filename);

%% read the tsv file

file = fopen(filename,'r');
headers = regexp(fgetl(file),'\t','split');

data = cell(0,length(headers));
while ~feof(file)
    data(end+1,:) = regexp(fgetl(file),'\t','split');
end

%% format the columns

t = table;
for i=1:length(headers)
    fprintf('\tReading %-15s', headers{i})
    temp = data(:,i);
    
    temp_num = cellfun2(@str2num,temp);
    if all(cellfun(@isempty,temp_num) | strcmpi(temp,'nan'))
        % string field
        if length(unique(temp))<.5*length(temp)
            % considered to be categorical
            fprintf('\t--> transformed in categorical\n')
            temp = categorical(temp);
        else
            fprintf('\t--> string field\n')
        end
    else
        if length(unique(cell2mat(temp_num)))==2 && all(unique(cell2mat(temp_num))==[0;1])
            fprintf('\t--> transformed in boolean\n')
            temp = cellfun(@logical,temp_num);
        else
            fprintf('\t--> transformed in numerical\n')
            temp = cell2mat(temp_num);
        end
    end
    t = [t table(temp, 'Variablenames', {headers{i}})];
end






