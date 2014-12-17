function t_data = Import_PlatesCellCountData(filename, plateinfo, varargin)
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
%           - Nuclei - Number Of Objects; can be changed to another field with
%               option 'NobjField', e.g. 
%               Import...(..., 'NobjField', 'Nuclei Selected - Number of Objects Nuclei Selected')
%           note: this allows to get multiple fields from the raw datafile;
%               the first field will be converted to Cellcount, other field
%               names will not be changed (except being made compatible).
%               For example: Import...(..., 'NobjField', ...
%                   {'Nuclei Selected - Number of Objects', ...
%                    'Nuclei - Number Of Objects'} )
%
%   Alternatively, filename can be a Nx[1-2] cell array with a filename and
%   corresponding barcode. Files will be read in order and if it is a Nx2
%   array, barcode will be added while reading. Each file should have
%   fields Well, NumberOfAnalyzedFields, Nuclei_NumberOfObjects, and a
%   field Result/Barcode if a Nx1 array is provided.
%
%   plateinfo should be a table (or the name of a tsv file) with headers:
%           - Barcode
%           - CellLine
%           - TreatmentFile (refer to a table or .hpdd/.mat with with the
%               plate design; use  -  if untreated, 1 if tab-separated file)
%           - DesignNumber (replicate or treatment design, mandatory for
%               .mat or .hpdd treatment files; set to 1 for other files)
%           - Replicate (optional - to differentiate to plates with the same
%               treatment file, and design number if .mat/.hpdd)
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
fprintf('Import Cell count data from Columbus files:\n');

p = inputParser;
addParameter(p,'NobjField',{'Nuclei_NumberOfObjects'},@(x) isstr(x) || iscellstr(x));
parse(p,varargin{:})
NobjField = p.Results.NobjField;
if ischar(NobjField), NobjField = {NobjField}; end
try
    NobjField = matlab.internal.tableUtils.makeValidName(NobjField,'silent');
catch
    NobjField = ReducName( ReplaceName(NobjField,['-/\:?!,.' 250:1e3], '_'), ' ');
end

%% check for proper inputs
if ischar(filename)
    assert(exist(filename,'file')>0, 'File %s is missing!', filename)
else
    assert(iscellstr(filename), 'Wrong type of input')
    for i=1:size(filename,1)
        assert(exist(filename{i,1},'file')>0, 'File %s is missing!', filename{i,1})
    end
end

if ischar(plateinfo)
    assert(exist(plateinfo,'file')>0, 'Barcode file %s is missing!', plateinfo)
    t_plateinfo = tsv2table(plateinfo);
elseif istable(plateinfo)
    t_plateinfo = plateinfo;
else
    error('Wrong argument for plateinfo')
end

CheckPlateInfo(t_plateinfo)
if iscell(t_plateinfo.DesignNumber)
    t_plateinfo.DesignNumber = cellstr2mat(t_plateinfo.DesignNumber);
end
%% load the cell count data

if ischar(filename)
    % case of a single file
    fprintf('\nImporting from: %s\n', filename)
    t_raw = tsv2table(filename);
else
    % case of a cell array (multiple files)
    t_raw = table;
    fprintf('\nImporting from: files\n')
    for i=1:size(filename,1)
        temp = tsv2table(filename{i,1});
        fprintf('\t%s\n', filename{i,1});
        
        if size(filename,2)==2
            % replace/ add the barcode
            if any(isvariable(temp, {'Result' 'Barcode'}))
                warnprintf(['Overwriting the Result/Barcode found in ' ...
                    'the individual file with input filename(:,2)'])
                if isvariable(temp, 'Barcode')
                    temp(:,'Barcode') =[];
                else
                    temp(:,'Result') =[];
                end
            end
            temp = [table(repmat(filename(i,2), height(temp),1), ...
                'VariableName', {'Barcode'}) temp];
        end
        
        t_raw = [t_raw; temp];
        
    end
end


% specific case of the output of Columbus: Results split into Barcode and
% date
if isvariable(t_raw, 'Result')
    t_raw = splitBarcodeDate(t_raw);
end

% correction of some variable names
CorrectionVarNames = {
    'WellName' 'Well'};
for i=1:size(CorrectionVarNames,1)
    if isvariable(t_raw, CorrectionVarNames{i,1})
        t_raw.Properties.VariableNames{CorrectionVarNames{i,1}} = ...
            CorrectionVarNames{i,2};
    end
end

% report error if missing field for cell count
if ~ismember(NobjField{1}, varnames(t_raw))
    warnprintf('Missing the field %s used for cell count', NobjField{1})
    disp('Available fields:')
    disp(varnames(t_raw)')
    error('Use optional input ''NobjField'' for specifying fields')
end
assert(all(ismember(NobjField, varnames(t_raw))), ...
    'Not all fields in ''NobjField'' are present in the file')
    

% check the number of fields
if length(unique(t_raw.NumberOfAnalyzedFields))>1
    Nref = median(t_raw.NumberOfAnalyzedFields);
    warnprintf('%i wells with missing fields', sum(t_raw.NumberOfAnalyzedFields<Nref))
    FieldCorrected = NobjField(strfindcell(NobjField,'Number')>0);
    warnprintf('Correcting for field number:\n\t - %s', strjoin(FieldCorrected,'\n\t - '))    
    for i = find(t_raw.NumberOfAnalyzedFields~=Nref)'
        for j=1:length(FieldCorrected)
            t_raw.(FieldCorrected{j})(i) = t_raw.(FieldCorrected{j})(i)*...
                (Nref/t_raw.NumberOfAnalyzedFields(i));
        end
    end
end


% correct well labels from \w%i to \w%02i
well_idx = cellfun(@length,t_raw.Well)==2;
t_raw.Well(well_idx) = strcat(cellcell2cellstr(regexp(t_raw.Well(well_idx),'^(\w)\d$','tokens')), ...
    '0',cellcell2cellstr(regexp(t_raw.Well(well_idx),'^\w(\d)$','tokens')));
%%
plate_barcodes = unique(t_raw.Barcode);

warnassert(length(plate_barcodes)>=height(t_plateinfo), ...
    'Found %i plates, but expected %i plates; missing %s', ...
    length(plate_barcodes), height(t_plateinfo), ...
    strjoin(setdiff(t_plateinfo.Barcode,plate_barcodes)',' ') );

otherVariables = setdiff(varnames(t_plateinfo), {'Time' 'CellLine' 'Barcode' ...
    'TreatmentFile' 'DesignNumber' 'ExpNumber'}, 'stable');
if ismember('Untrt', otherVariables)
    if ~all(strcmp(t_plateinfo.TreatmentFile(t_plateinfo.Untrt>0),'-')) || ...
            ~all(t_plateinfo.Untrt(strcmp(t_plateinfo.TreatmentFile,'-'))>0)
        warnprintf('Discrepency between TreatmentFile==''-'' and Untrt; Overwriting TreatmentFile')
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
    [~,~,ext] = fileparts(t_plateinfo.TreatmentFile{iBC});
    
    if (strcmp(ext,'.mat') || strcmp(ext,'.hpdd'))
        assert(isvariable(t_plateinfo,'DesignNumber'), ...
            'Treatment files .mat or .hpdd needs a DesignNumber column')
        DesignNumber(idx) = t_plateinfo.DesignNumber(iBC);
    elseif ~isempty(ext) % assume a tab-separated file
        if isvariable(t_plateinfo,'DesignNumber') && ...
                (t_plateinfo.DesignNumber(iBC)~=1 || ...
                ~isnan(t_plateinfo.DesignNumber(iBC)))
            warnprintf('for tab-separated file, DesignNumber is forced to be 1')
        end
        DesignNumber(idx) = 1;
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
Untrt = cellfun(@(x) strcmp(x,'-') || isempty(x), TreatmentFile) ;
assert(~any(DesignNumber==0 & ~Untrt), 'Some wells are not ''Untrt'' and don''t have a DesignNumber')
warnassert(all(Time(~Untrt)>0), 'Some treated wells don''t have a Time')

% compile the finale table
t_data = [table(Barcode, CellLine, TreatmentFile, DesignNumber, Untrt, Time) ...
    t_raw(Usedidx, intersect([{'Well'} NobjField {'Date'}], varnames(t_raw), 'stable'))];
if ~isempty(otherVariables)
    fprintf(['\tAdded variable(s): ''' cellstr2str(otherVariables, ''', ''') '''\n'])
    eval(['t_data = [t_data table(' cellstr2str(otherVariables, ',') ')];'])
end
t_data.Properties.VariableNames{NobjField{1}} = 'Cellcount';
t_data = TableToCategorical(t_data);

fprintf('\n')
end


function t_raw = splitBarcodeDate(t_raw)
Code_date = regexp(t_raw.Result,' > ','split');
Code_date = vertcat(Code_date{:});
t_raw = [cell2table(Code_date,'VariableName',{'Barcode' 'Date'}) t_raw(:,2:end)];
end

function CheckPlateInfo(t_plateinfo)
Infovars = {'Barcode' 'Time' 'CellLine' 'TreatmentFile'};
for i=1:length(Infovars)
    assert(isvariable(t_plateinfo, Infovars{i}), 'Missing columns %s in the plate info',...
        Infovars{i})
end
for i=1:height(t_plateinfo)
    tf = t_plateinfo.TreatmentFile{i};
    [~,n,e] = fileparts(tf);
    assert(ismember(e, {'.mat' '.hpdd' '.tsv' '.txt'}) | strcmp(n,'-'), ...
        'TreatmentFile should be a .mat, .hpdd, .txt, .tsv file of ''-'' for untreated')
end
end



