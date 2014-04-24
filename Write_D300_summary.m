function Write_D300_summary(filename)

%%

[~,sheets] = xlsfinfo(filename);
trt_sheets = sheets(strfindcell(sheets,'_total')>0);

%%

allDrugs = cell(0,1);
stock = cell(0,1);
volumes = zeros(0,length(trt_sheets));


for iS = 1:length(trt_sheets)

[~,~,data] = xlsread(filename,trt_sheets{iS});

%
rowidx = find(strcmp(data(:,2),'Stock (mM)='));
%%
for iR=1:length(rowidx)
    drug = data{rowidx(iR),1};
    idx = strcmp(allDrugs,drug);
    if any(idx)
        assert(stock{idx}==data{rowidx(iR),3})
        volumes(idx,iS) = data{rowidx(iR),5};
    else
        allDrugs = [allDrugs; drug];
        stock = [stock; data(rowidx(iR),3)];
        volumes(length(allDrugs),iS) = data{rowidx(iR),5};
    end
end

%
end 
allDrugs = regexp(allDrugs,' ', 'split');
allDrugs = vertcat(allDrugs{:});

%%


[~,sheets] = xlsfinfo(filename);
if ismember('Trt_summary', sheets)
    [~,~,dumb] = xlsread(filename, 'Trt_summary');
    xlswrite(filename, cell(size(dumb)), 'Trt_summary');
end

xlswrite(filename, [{'DrugName' 'HMSL id' 'Stock (mM)'} strcat(trt_sheets, '_(ul)') {'round-up Sum (ul)'};
    allDrugs stock num2cell(round(volumes)/1e3) num2cell(ceil((sum(volumes,2)+.05)/1e2)/10)], 'Trt_summary')



