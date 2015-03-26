function SingleCell = Import_SingleCell(t_processed, folder)
% SingleCell = Import_SingleCell(t_processed, folder)
%   t_processed     table outputed from Merge_CellCountData; need the
%                       columns 'Barcode' 'Well' and will condiser the data
%                       columns that ends with '_MeanPerWell'
%
%   folder          folder in which the single cell data are stored
%
%   SingleCell      structure with the data fields for each plate/well
%                   

%%
meanfields = varnames(t_processed)';
meanfields = meanfields(regexpcell(meanfields, '_MeanPerWell')>0);
meanfields = meanfields(regexpcell(meanfields, 'Ctrl_')==0);
meanfields = meanfields(regexpcell(meanfields, 'WholeImage')==0);

meanfields = regexpcelltokens(meanfields, '(\w*)_MeanPerWell',1);

SingleCell = struct;
for i=1:length(meanfields)
    SingleCell.(meanfields{i}) = [];
end
SingleCell.Barcode = ''; SingleCell.Well = '';
SingleCell = repmat(SingleCell, height(t_processed),1);
%%

for iPW = 1:height(t_processed)
    
    
    folderlist = dir(folder);
    allfolder = {folderlist([folderlist.isdir]==1).name};
    
    Plate = char(t_processed.Barcode(iPW));
    Well = char(t_processed.Well(iPW));
    SingleCell(iPW).Barcode = Plate;
    SingleCell(iPW).Well = Well;
    
    if Well(2)=='0', Well = Well([1 3]); end
    
    subfolder = [folder filesep allfolder{regexpcell(allfolder, Plate)>0}];
    
    files = dir(subfolder);
    files = {files([files.isdir]==0).name};
    files = files(regexpcell(files, 'Whole Image')==0);
    file = files(regexpcell(files, sprintf('result.%s\\[', Well))>0);
    assert(length(file)==1);
    
    %%
    t_ss = tsv2table([subfolder filesep file{1}]);
    
    fprintf('\nReading Plate %s, well %s', Plate, Well)
    
    for i=1:length(meanfields)
        SingleCell(iPW).(meanfields{i}) = t_ss.(meanfields{i});
    end
    
end


