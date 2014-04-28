function Write_D300_summary(filename)
%%

[~,sheets] = xlsfinfo(filename);

%%
sheets_total = sheets(strfindcell(sheets,'_total')>0);

Drugs = {};
Stock = [];
Volume = [];

output = ['Drug Name' 'HMSL id' 'Stock (mM)' ...
    strcat(sheets_total, ' (ul)') 'Total volume (ul)'];

for iS = 1:length(sheets_total)
    
    [~,~,data] = xlsread(filename, sheets_total{iS});
    
    idx = find(cellfun(@isstr, data(:,2)));
    idx = idx(strfindcell(data(idx,2),'Stock (mM)')==1);

    for i=1:length(idx)
        t_Drug = data{idx(i),1};
        Didx = find(strcmp(Drugs, t_Drug));
    
        if isempty(Didx)
            Drugs{end+1} = t_Drug;
            Stock(end+1) = data{idx(i),3};
            Volume(end+1,iS) = data{idx(i),5};
            if strfind(t_Drug, 'HMSL')
                t_Drug = regexp(t_Drug,'HMSL','split');
                output(end+1,1:3) = {t_Drug{1} ['HMSL' t_Drug{2}] num2str(data{idx(i),3})};
            else
                output(end+1,1:3) = {t_Drug '_' num2str(data{idx(i),3})};
            end
            output{end,iS+3} = num2str(data{idx(i),5}/1e3,'%.2f');
        else
            assert(Stock(Didx)==data{idx(i),3}, 'not the same stock conc for %s', ...
                Drugs{Didx});
            Volume(Didx,iS) = data{idx(i),5};
            output{Didx+1,iS+3} = num2str(data{idx(i),5}/1e3,'%.2f');
        end
    end
end

for i=1:length(Drugs)
    output{i+1,end} = num2str(ceil(sum(Volume(i,:)/100))/10,'%.1f');
end
            

    
%%
sheetname = 'Summary_drugs';

if ismember(sheetname, sheets)
    fprintf('Replacing sheet %s in file %s\n', sheetname, filename)
    [~,~,dumb] = xlsread(filename, sheetname);
    xlswrite(filename, cell(size(dumb)), sheetname);
end


xlswrite(filename, output, sheetname)