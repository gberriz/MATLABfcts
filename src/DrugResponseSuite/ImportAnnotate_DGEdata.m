
function [mRNA_values, mRNA_labels] = ImportAnnotate_DGEdata(filename, plateinfo, Barcode, folder)

if ~strcmp(filename((end-3):end), '.dat')
    filename = [filename '.unq.refseq.umi.dat'];
end

p = struct('TimeCourse', false, 'Cellcount', 'none');
t_plateinfo = ImportCheckPlateInfo(plateinfo, p);

rawRNA = importdata(filename);


%%

Wells = CorrectWellsTo0xY(regexpcelltokens(rawRNA.textdata(1,2:end), '_([A-Q][0-9]*)'))';
t_labels = table(repmat({Barcode}, length(Wells),1), Wells, ...
    'variablenames', {'Barcode' 'Well'});

t_labels = AddPlateInfo_RawData(t_labels, t_plateinfo, '', p);

t_labels = Annotate_CellCountData(t_labels, folder);

mRNA_names = table(categorical(rawRNA.textdata(2:end,1)),'VariableNames',{'gene'});
mRNA_values = rawRNA.data;

mRNA_labels = {mRNA_names t_labels};
