function t_data = Import_PlatesCellCountData(filename, plateinfo, varargin)
% t_data = Import_PlatesCellCountData(filename, plateinfo, varargin)
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
%           - Time (in hours); 
%                   note: Time is ignored if the option 'TimeCourse' is
%                   true; in that case, the timestamp in the Columbus file
%                   will be used to tag the time of the sample. 
%           - any Additional field with relevant plate properties to pass to the
%               data (collagen, ...)
%           note: 'ExpNumber' is ignored as it serves as an internal
%               control for the D300
%
%   varargin:   - 'NobjField'   [Nuclei - Number Of Objects]
%               - 'Cellcount'   [Nuclei - Number Of Objects]
%               - 'TimeCourse'  [false]
%               - 'T0shift'     [0.25 h] Different in hours between the
%                                   first plate imaged and the treatment
%               - 'T0date'      Input alternative to T0shift for timecourse:
%                                   date and time of the treatment
%
%   t_data is a table with each well annotated accorting to the barcode.
%   The column 'Untrt' is evaluated and the data are corrected for the
%   number of fields (stored in column Cellcount).
%       the column Cellcount can be a function of the input columns (by
%       default Cellcount is 'Nuclei_NumberOfObjects' or the first column
%       in the field 'NobjField').
%
%   Example:
%   --------
%
%   t_data = Import_PlatesCellCountData('Results_20150320.txt'], ...
%     'PlateIDs_20150320.txt', 'NobjField', ...
%     {'Nuclei-Hoechst - Number of Objects' 'Nuclei-LDRpos - Number of Objects'},...
%     'Cellcount', @(x) x(:,1)-x(:,2));
%

fprintf('Import Cell count data from Columbus files:\n');

p = inputParser;
addParameter(p,'NobjField',{'Nuclei_NumberOfObjects'},@(x) ischar(x) || iscellstr(x));
addParameter(p,'Cellcount', [], @(x) isa(x,'function_handle'));
addParameter(p,'TimeCourse', false, @islogical);
addParameter(p,'T0shift', 1/4, @isscalar);
addParameter(p,'T0date', [], @isvector);
parse(p,varargin{:})
p = p.Results;

NobjField = p.NobjField;
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

t_plateinfo = ImportCheckPlateInfo(plateinfo, p);

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
if length(NobjField)>1 && ~strcmp(NobjField{2},'all')
    if ~all(ismember(NobjField, varnames(t_raw)))
        warnprintf('Not all fields in ''NobjField'' are present in the file')
        disp('Available fields:')
        disp(varnames(t_raw)')
        error('Missing fields: %s', strjoin(NobjField(~ismember(NobjField, varnames(t_raw)))))
    end
elseif length(NobjField)==2 && strcmp(NobjField{2},'all')
    NobjField = [NobjField(1) setdiff(varnames(t_raw), [NobjField(1) ...
        {'Well' 'Barcode' 'Result' 'Date' 'Row' 'Column' ...
        'NumberOfAnalyzedFields' 'URL'}])];
end

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

t_raw.Well = CorrectWellsTo0xY(t_raw.Well);
% % correct well labels from \w%i to \w%02i
% well_idx = cellfun(@length,t_raw.Well)==2;
% t_raw.Well(well_idx) = strcat(cellcell2cellstr(regexp(t_raw.Well(well_idx),'^(\w)\d$','tokens')), ...
%     '0',cellcell2cellstr(regexp(t_raw.Well(well_idx),'^\w(\d)$','tokens')));
%%
t_data = AddPlateInfo_RawData(t_raw, t_plateinfo, NobjField, p);

% plate_barcodes = unique(t_raw.Barcode);
% 
% warnassert(length(plate_barcodes)>=height(t_plateinfo), ...
%     'Found %i plates, but expected %i plates; missing %s', ...
%     length(plate_barcodes), height(t_plateinfo), ...
%     strjoin(setdiff(t_plateinfo.Barcode,plate_barcodes)',' ') );
% 
% otherVariables = setdiff(varnames(t_plateinfo), {'Time' 'CellLine' 'Barcode' ...
%     'TreatmentFile' 'DesignNumber' 'ExpNumber'}, 'stable');
% if ismember('Untrt', otherVariables)
%     if ~all(strcmp(t_plateinfo.TreatmentFile(t_plateinfo.Untrt>0),'-')) || ...
%             ~all(t_plateinfo.Untrt(strcmp(t_plateinfo.TreatmentFile,'-'))>0)
%         warnprintf('Discrepency between TreatmentFile==''-'' and Untrt; Overwriting TreatmentFile')
%         t_plateinfo.TreatmentFile(t_plateinfo.Untrt>0) = {'-'};
%     end
%     t_plateinfo.Untrt = [];
%     otherVariables = setdiff(otherVariables, 'Untrt', 'stable');
% end
% 
% CellLine = cell(height(t_raw),1);
% Barcode = cell(height(t_raw),1);
% TreatmentFile = cell(height(t_raw),1);
% DesignNumber = zeros(height(t_raw),1);
% Time = zeros(height(t_raw),1);
% Untrt = false(height(t_raw),1);
% 
% 
% cnt = 0;
% % join the barcode table to t_raw with a few controls
% for iBC = 1:height(t_plateinfo)
%     
%     idx = find(strcmp(t_raw.Barcode, t_plateinfo.Barcode{iBC}));
%     if isempty(idx)
%         warnprintf('No result found for plate %s', t_plateinfo.Barcode{iBC})
%     end
%     
%     CellLine(idx) = t_plateinfo.CellLine(iBC);
%     Barcode(idx) = t_plateinfo.Barcode(iBC);
%     TreatmentFile(idx) = t_plateinfo.TreatmentFile(iBC);
%     [~,~,ext] = fileparts(t_plateinfo.TreatmentFile{iBC});
%     
%     if (strcmp(ext,'.mat') || strcmp(ext,'.hpdd'))
%         assert(isvariable(t_plateinfo,'DesignNumber'), ...
%             'Treatment files .mat or .hpdd needs a DesignNumber column')
%         DesignNumber(idx) = t_plateinfo.DesignNumber(iBC);
%     elseif ~isempty(ext) % assume a tab-separated file
%         if isvariable(t_plateinfo,'DesignNumber') && ...
%                 (t_plateinfo.DesignNumber(iBC)~=1 || ...
%                 ~isnan(t_plateinfo.DesignNumber(iBC)))
%             warnprintf('for tab-separated file, DesignNumber is forced to be 1')
%         end
%         DesignNumber(idx) = 1;
%     end
%     Untrt(idx) = strcmp(t_plateinfo.TreatmentFile(iBC),'-');
%     if p.TimeCourse
%         % rounded to the 1/100 of hour, 1st value is p.T0shift (default ~15 minutes).
%         if ~isempty(p.T0date)
%             Time(idx) = .01*round(100*( 24*(t_raw.Date(idx)-datenum(p.T0date) )));
%         else
%             Time(idx) = .01*round(100*( 24*(t_raw.Date(idx)-min(t_raw.Date)) + p.T0shift ));
%         end
%     else
%         Time(idx) = t_plateinfo.Time(iBC);
%     end
%     
%     % parse the additional plate information from the barcode file
%     for i = 1:length(otherVariables)
%         eval([otherVariables{i} '(idx) = t_plateinfo.' otherVariables{i} '(iBC);'])
%     end
%     
%     cnt = cnt+1;
% end
% 
% % format properly the additional plate information from the plateinfo file
% for i = 1:length(otherVariables)
%     eval([otherVariables{i} ' = ToColumn(' otherVariables{i} ');'])
%     eval(['if length(' otherVariables{i} ')<height(t_raw), ' ...
%         otherVariables{i} '(height(t_raw)+1) = ' otherVariables{i} '(1);' ...
%         otherVariables{i} '(height(t_raw)+1) = []; end'])
% end
% 
% % broadcast the properties
% if cnt<length(plate_barcodes)
%     warning('some entries in the result file are unused!')
%     Usedidx = ~cell2mat(cellfun2(@isempty,CellLine));
%     TreatmentFile = TreatmentFile(Usedidx);
%     CellLine = CellLine(Usedidx);
%     DesignNumber = DesignNumber(Usedidx);
%     Time = Time(Usedidx);
%     Barcode = Barcode(Usedidx);
%     for i = 1:length(otherVariables)
%         eval([otherVariables{i} ' = ' otherVariables{i} '(Usedidx);'])
%     end
% else
%     Usedidx = 1:height(t_raw);
% end
% 
% % Evaluate the untreated plates
% Untrt = cellfun(@(x) strcmp(x,'-') || isempty(x), TreatmentFile) ;
% assert(~any(DesignNumber==0 & ~Untrt), 'Some wells are not ''Untrt'' and don''t have a DesignNumber')
% warnassert(all(Time(~Untrt)>0), 'Some treated/perturbed wells have Time=0')
% 
% % compile the final table
% t_data = [table(Barcode, CellLine, TreatmentFile, DesignNumber, Untrt, Time) ...
%     t_raw(Usedidx, intersect([{'Well'} NobjField {'Date'}], varnames(t_raw), 'stable'))];
% if ~isempty(otherVariables)
%     fprintf(['\tAdded variable(s): ''' cellstr2str(otherVariables, ''', ''') '''\n'])
%     eval(['t_data = [t_data table(' cellstr2str(otherVariables, ',') ')];'])
% end
% if isempty(p.Cellcount)
%     t_data.Properties.VariableNames{NobjField{1}} = 'Cellcount';
% else
%     temp = table2array(t_data(:,NobjField));
%     temp = p.Cellcount(temp);
%     t_data.Cellcount = temp;
% end
% t_data = TableToCategorical(t_data);

fprintf('\n')
end


function t_raw = splitBarcodeDate(t_raw)
Code_date = regexp(t_raw.Result,' > ','split');
Code_date = vertcat(Code_date{:});
Code_date(:,2) = regexpcellsplit(Code_date(:,2), '(',1);
Code_date(:,2) = cellfun2(@datenum, Code_date(:,2));
t_raw = [cell2table(Code_date,'VariableName',{'Barcode' 'Date'}) t_raw(:,2:end)];
end
