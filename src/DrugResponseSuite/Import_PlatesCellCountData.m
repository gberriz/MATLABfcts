function t_data = Import_PlatesCellCountData(filename, plateinfo)
% t_data = Import_PlatesCellCountData(filename, plateinfo)
%
%   merge the information in the barcode (either file or table) with the
%   data from the output of Columbus or from other tsv file with proper
%   columns.
%
%   filename is a tsv file which should contain at least the headers
%   (default output from Columbus):
%           - Result (containing the barcodes and the timestamp)
%               -> or two columns Barcode and Date
%           - Well
%           - NumberOfAnalyzedFields
%           - Nuclei_NumberOfObjects
%           note: cannot currently handle other columns
%
%   plateinfo should be a table (or the name of a tsv file) with headers:
%           - Barcode
%           - CellLine
%           - TreatmentFile (refer to a table or .hpdd/.mat with with the
%               plate design; use  -  if untreated)
%           - DesignNumber (replicate or treatment design, mandatory for
%               .mat or .hpdd treatment files)
%           - Time (in hours)
%           - any Additional field with relevant plate properties to pass to the
%               data (collagen, ...)
%           note: 'ExpNumber' is ignored as it serves as an internal
%               control for the D300
%
%   t_data is a table with each well annotated accorting to the barcode.
%   The column 'Untrt' is evaluated and the data are corrected for the
%   number of fields (stored in column Cellcount).
%
%


%% check for proper inputs
assert(exist(filename,'file')>0, 'File %s is missing!', filename)

if ischar(plateinfo)
    assert(exist(plateinfo,'file'), 'Barcode file %s is missing!', plateinfo)
    t_plateinfo = tsv2table(plateinfo);
elseif istable(plateinfo)
    t_plateinfo = plateinfo;
else
    error('Wrong argument for plateinfo')
end

fprintf('\nImporting from: %s\n', filename)
%% load the cell count data
t_raw = shutup(@() tsv2table(filename));

% specific case of the output of Columbus: Results split into Barcode and
% date
if isvariable(t_raw, 'Result')
    t_raw = splitBarcodeDate(t_raw);
end

% check the number of fields
if length(unique(t_raw.NumberOfAnalyzedFields))>1
    Nref = median(t_raw.NumberOfAnalyzedFields);
    fprintf('\tWarning: %i wells with missing fields\n', ...
        sum(t_raw.NumberOfAnalyzedFields<Nref))
    for i = find(t_raw.NumberOfAnalyzedFields~=Nref)'
        t_raw.Nuclei_NumberOfObjects(i) = t_raw.Nuclei_NumberOfObjects(i)*...
            (Nref/t_raw.NumberOfAnalyzedFields(i));
    end
end


% correct well labels from \w%i to \w%02i
well_idx = cellfun(@length,t_raw.Well)==2;
t_raw.Well(well_idx) = strcat(cellcell2cellstr(regexp(t_raw.Well(well_idx),'^(\w)\d$','tokens')), ...
    '0',cellcell2cellstr(regexp(t_raw.Well(well_idx),'^\w(\d)$','tokens')));
%%
plate_barcodes = unique(t_raw.Barcode);

assert(length(plate_barcodes)>=height(t_plateinfo), ...
    'Found %i plates, but expected %i plates; missing %s', ...
    length(plate_barcodes), height(t_plateinfo), ...
    strjoin(setdiff(t_plateinfo.Barcode,plate_barcodes)',' ') );

otherVariables = setdiff(varnames(t_plateinfo), {'Time' 'CellLine' 'Barcode' ...
    'TreatmentFile' 'DesignNumber' 'ExpNumber'}, 'stable');
if ismember('Untrt', otherVariables)
    if ~all(strcmp(t_plateinfo.TreatmentFile(t_plateinfo.Untrt>0),'-')) || ...
            ~all(t_plateinfo.Untrt(strcmp(t_plateinfo.TreatmentFile,'-'))>0)
        warning('Discrepency between TreatmentFile==''-'' and Untrt; Overwriting TreatmentFile')
        t_plateinfo.TreatmentFile(t_plateinfo.Untrt>0) = {'-'};
    end
    t_plateinfo.Untrt = [];
    otherVariables = setdiff(otherVariables, 'Untrt', 'stable');
end

CellLine = cell(height(t_raw),1);
Barcode = cell(height(t_raw),1);
TreatmentFile = cell(height(t_raw),1);
DesignNumber = zeros(height(t_raw),1);
Time = zeros(height(t_raw),1);
Untrt = false(height(t_raw),1);


cnt = 0;
% join the barcode table to t_raw with a few controls
for iBC = 1:height(t_plateinfo)
    
    idx = find(strcmp(t_raw.Barcode, t_plateinfo.Barcode{iBC}));
    
    CellLine(idx) = t_plateinfo.CellLine(iBC);
    Barcode(idx) = t_plateinfo.Barcode(iBC);
    TreatmentFile(idx) = t_plateinfo.TreatmentFile(iBC);
    if ~isempty(strfind(t_plateinfo.TreatmentFile{iBC}, '.mat')) || ...
            isvariable(t_plateinfo,'DesignNumber')
        assert(isvariable(t_plateinfo,'DesignNumber'), ...
            'Treatment files .mat needs a DesignNumber column')
        DesignNumber(idx) = t_plateinfo.DesignNumber(iBC);
    else
        DesignNumber(idx) = NaN;
    end
    Untrt(idx) = strcmp(t_plateinfo.TreatmentFile(iBC),'-');
    Time(idx) = t_plateinfo.Time(iBC);
    
    % parse the additional plate information from the barcode file
    for i = 1:length(otherVariables)
        eval([otherVariables{i} '(idx) = t_plateinfo.' otherVariables{i} '(iBC);'])
    end
    
    cnt = cnt+1;    
end

% format properly the additional plate information from the plateinfo file
for i = 1:length(otherVariables)
    eval([otherVariables{i} ' = ToColumn(' otherVariables{i} ');'])
    eval(['if length(' otherVariables{i} ')<height(t_raw), ' ...
        otherVariables{i} '(height(t_raw)+1) = ' otherVariables{i} '(1);' ...
        otherVariables{i} '(height(t_raw)+1) = []; end'])
end

% broadcast the properties
assert(cnt<=length(plate_barcodes))
if cnt<length(plate_barcodes)
    warning('some entries in the result file are unused!')
    Usedidx = ~cell2mat(cellfun2(@isempty,CellLine));
    TreatmentFile = TreatmentFile(Usedidx);
    CellLine = CellLine(Usedidx);
    DesignNumber = DesignNumber(Usedidx);
    Time = Time(Usedidx);
    for i = 1:length(otherVariables)
        eval([otherVariables{i} ' = ' otherVariables{i} '(Usedidx);'])
    end
else
    Usedidx = 1:height(t_raw);
end

% Evaluate the untreated plates
Untrt = cellfun(@(x) strcmp(x,'-') || isempty(x), TreatmentFile);
assert(~any(DesignNumber==0 & ~Untrt), 'Some wells are not ''Untrt'' and don''t have a DesignNumber')
assert(all(Time(~Untrt)>0), 'Some treated wells don''t have a Time')

% compile the finale table
t_data = [table(CellLine,TreatmentFile,DesignNumber,Untrt,Time) ...
    t_raw(Usedidx, {'Well' 'Row' 'Column' 'Nuclei_NumberOfObjects' 'Date'})];
if ~isempty(otherVariables)
    fprintf(['\tAdded variable(s) ''' cellstr2str(otherVariables, ''', ''') '''\n'])
    eval(['t_data = [t_data table(' cellstr2str(otherVariables, ',') ')];'])
end
t_data.Properties.VariableNames{'Nuclei_NumberOfObjects'} = 'Cellcount';
t_data = TableToCategorical(t_data,[1 2 6]);

fprintf('\n')
end


function t_raw = splitBarcodeDate(t_raw)
Code_date = regexp(t_raw.Result,' > ','split');
Code_date = vertcat(Code_date{:});
t_raw = [cell2table(Code_date,'VariableName',{'Barcode' 'Date'}) t_raw(:,2:end)];
end