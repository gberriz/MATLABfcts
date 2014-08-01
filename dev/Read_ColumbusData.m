function t_data = Read_ColumbusData(filename, barcode_file)
% t_data = Read_ColumbusData(filename, barcode_file)
%
%   import the results of plates scanned with Columbus
%
%   filename is the output of Columbus which should contain at least the headers:
%           Result (matching the barcodes)
%           Well
%           NumberOfAnalyzedFields
%           Nuclei_NumberOfObjects
%           URL (for checking the PlateID)
%
%   barcode_file should be a table with experimental design with headers:
%           Barcode
%           CellLine
%           Treatmentfile (refer to a table or .mat with with the
%               treatments on the plate - plate design) use   -   if untreated
%           Replicate:  replicate (treatment)
%           time (in hours)
%

t_raw = readtsv(filename);


if length(unique(t_raw.NumberOfAnalyzedFields))>1    
    Nref = median(t_raw.NumberOfAnalyzedFields);
    warning('wells with missing fields: %i', ...
        unique(t_raw.NumberOfAnalyzedFields))
    for i = find(t_raw.NumberOfAnalyzedFields~=Nref)'
        t_raw.Nuclei_NumberOfObjects(i) = t_raw.Nuclei_NumberOfObjects(i)*...
            (Nref/t_raw.NumberOfAnalyzedFields(i));
    end
end

Code_date = regexp(t_raw.Result,' > ','split');
Code_date = vertcat(Code_date{:});
t_raw = [cell2table(Code_date,'VariableName',{'Barcode' 'date'}) t_raw(:,2:end)];

%%

t_barcode = readtsv(barcode_file); 

assert(length(unique(t_raw.Barcode))>=height(t_barcode), ...
    'Found %i plates, but expected %i plates; missing %s', ...
    length(unique(t_raw.Barcode)), height(t_barcode), ...
    strjoin(setdiff(t_barcode.Barcode,unique(t_raw.Barcode))',' ') );

url = regexp(t_raw.URL,'/','split');
PlateID = vertcat(url{:}); 

PlateID = cellfun(@str2num,PlateID(:,end-1));

CellLine = cell(height(t_raw),1);
Barcodes = cell(height(t_raw),1);
Treatmentfile = cell(height(t_raw),1);
Replicate = zeros(height(t_raw),1);
Time = zeros(height(t_raw),1);

cnt = 0;
for iBC = 1:height(t_barcode)
    
    idx = find(strcmp(t_raw.Barcode, t_barcode.Barcode{iBC}));
    assert(length(unique(PlateID(idx)))==1, '%i plate IDs found', ...
        length(unique(PlateID(idx))))
    assert(isempty(intersect(unique(PlateID(idx)), unique(PlateID(~idx)))))
    
    CellLine(idx) = t_barcode.CellLine(iBC);
    Barcodes(idx) = t_barcode.Barcode(iBC);
    Treatmentfile(idx) = t_barcode.Treatmentfile(iBC);
    Replicate(idx) = t_barcode.Replicate(iBC);
        
    Time(idx) = t_barcode.Time(iBC);
    
    cnt = cnt+1;
    
end
assert(cnt<=length(unique(PlateID)))
if cnt<length(unique(PlateID))
    warning('some entries in the result file are unused!')
    Usedidx = ~cell2mat(cellfun2(@isempty,CellLine));
    Treatmentfile = Treatmentfile(Usedidx);
    CellLine = CellLine(Usedidx);
    Replicate = Replicate(Usedidx);
    Untrt = Untrt(Usedidx);
    Time = Time(Usedidx);    
else
    
    Usedidx = 1:height(t_raw);
end

Untrt = cellfun(@(x) strcmp(x,'-') || isempty(x), Treatmentfile);
assert(all(Replicate>0 | Untrt))
assert(all(Time(~Untrt)>0))

t_data = [table(CellLine,Treatmentfile,Replicate,Untrt,Time) t_raw(Usedidx, ...
    {'Well' 'Row' 'Column' 'Nuclei_NumberOfObjects' 'date'})];
t_data.Properties.VariableNames{'Nuclei_NumberOfObjects'} = 'Cellcount';
t_data = TableToCategorical(t_data,[1 2 6]);