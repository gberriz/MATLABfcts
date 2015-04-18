function SingleCell = Import_SingleCell(t_processed, folder, timecourse, filefilter)
% SingleCell = Import_SingleCell(t_processed, folder, timecourse)
%   t_processed     table outputed from Merge_CellCountData; need the
%                       columns 'Barcode' 'Well' and will condiser the data
%                       columns that ends with '_MeanPerWell'
%                       If column 'Date' exists, it will assume a timecourse 
%                       and look for corresponding file name (unless
%                       variable timecourse is false)
%
%   folder          folder in which the single cell data are stored
%
%   timecourse      time course experiment, multiple time points for the
%                       same well (optional, default=false)
%   
%   filefilter      regular experssion for filtering files in case of
%                       multiple output files for the same plate/well
%
%   SingleCell      structure with the data fields for each plate/well
%                   

%%
if ~exist('timecourse', 'var') || isempty(timecourse)
    timecourse = isvariable(t_processed,'Date');
end

meanfields = varnames(t_processed)';
meanfields = meanfields(regexpcell(meanfields, '_MeanPerWell')>0);
meanfields = meanfields(regexpcell(meanfields, 'Ctrl_')==0);
meanfields = meanfields(regexpcell(meanfields, 'WholeImage')==0);

meanfields = regexpcelltokens(meanfields, '(\w*)_MeanPerWell',1);

SingleCell = struct;
for i=1:length(meanfields)
    SingleCell.(meanfields{i}) = [];
end
SingleCell.Barcode = ''; SingleCell.Well = ''; SingleCell.Date = '';
SingleCell = repmat(SingleCell, height(t_processed),1);
%%

assert(exist(folder,'dir')>0, ['Folder: ' folder ' missing'])

t_out = t_processed;
output_cnt = 0;
for iPW = 1:height(t_processed)    
    
    folderlist = dir(folder);
    allfolder = {folderlist([folderlist.isdir]==1).name};
    
    Plate = char(t_processed.Barcode(iPW));
    Well = char(t_processed.Well(iPW));
    
    if Well(2)=='0', Well = Well([1 3]); end
    
    subfolder = [folder filesep allfolder{regexpcell(allfolder, Plate)>0}];
    
    files = dir(subfolder);
    files = {files([files.isdir]==0).name};
    files = files(regexpcell(files, 'Whole Image')==0);
    files = files(regexpcell(files, sprintf('result.%s\\[', Well))>0);
    
    if isvariable(t_processed,'Date') && timecourse
        files = files(regexpcell(files, ...
            datestr(t_processed.Date(iPW),'yyyy-mm-ddTHHMMSS'))==1);
    end
    
    if exist('filefilter','var')
        files = files(regexpcell(files, filefilter)>0);
    end
    
    if isempty(files)
        error(['Missing file for: ' Plate ' ' Well]);
    end
    
    if ~timecourse || isvariable(t_processed,'Date')
       assert(length(files)==1, ['Too many files: ' strjoin(files,' - ')])
    end
    
    for it = 1:length(files)
        date = regexp(files{it},'^([0-9\-]*)T','tokens');
        time = regexp(files{it}, '^[0-9\-]*T([0-9]*)\-', 'tokens');
        Date = datenum([date{1}{1} '-' time{1}{1}], 'yyyy-mm-dd-HHMMSS');
        
        fprintf('\nReading Plate %s, well %s, date %s', ...
            Plate, Well, [date{1}{1} '-' time{1}{1}])
        
        t_ss = tsv2table([subfolder filesep files{it}]);
        
        output_cnt = output_cnt+1;
        
        if isvariable(t_processed,'Date')
            t_out(output_cnt,:) = t_processed(iPW,:);
        else
            t_out(output_cnt,:) = [t_processed(iPW,:) table(Date)];
        end
        SingleCell(output_cnt).Barcode = Plate;
        SingleCell(output_cnt).Well = Well;
        SingleCell(output_cnt).Date = Date;
        for i=1:length(meanfields)
            SingleCell(output_cnt).(meanfields{i}) = t_ss.(meanfields{i});
        end
    end
    
end

fprintf('\n\n');

