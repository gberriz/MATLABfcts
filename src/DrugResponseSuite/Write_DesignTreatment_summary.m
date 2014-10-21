function Write_DesignTreatment_summary(filename, Designs, t_plateinfo)
% Write_D300_summary(filename, Designs, t_plateinfo)

% load the table is a file name was passed
if ischar(t_plateinfo)
    t_plateinfo = tsv2table(t_plateinfo);
else
    assert(istable(t_plateinfo), ['Barcodes should be a table or a tsv file' ...
        ' with columns Barcode, DesignNumber'])
end

%%

Drugs = get_DesignDrugs( Designs(setdiff(unique(t_plateinfo.DesignNumber),0)) );

output = ['Drug Name' 'HMSL id' 'Stock (mM)' ...
    strcat(t_plateinfo.Barcode)' 'D300 Total volume (ul)'];

Volume = zeros(length(Drugs), height(t_plateinfo));

for iPlte = 1:height(t_plateinfo)
    DesIdx = t_plateinfo.DesignNumber(iPlte);
    if DesIdx==0
        continue
    end
    
    WellVolume = Designs(DesIdx).well_volume;
    for iDr=1:length(Designs(DesIdx).Drugs)
        
        Drug = Designs(DesIdx).Drugs(iDr).DrugName;
        Didx = find(strcmp({Drugs.DrugName}, Drug));
    
        stock_conc = Drugs(Didx).stock_conc;
        
        Volume(Didx,iPlte) = ceil(100*sum(Designs(DesIdx).Drugs(iDr).layout(:)) ...
            *WellVolume/stock_conc)/100;        
    end
end

output(2:(1+length(Drugs)),1) = {Drugs.name}';
output(2:end,2) = {Drugs.HMSLid}';
output(2:end,3) = cellfun(@(x) {x/1e3}, {Drugs.stock_conc}');
output(2:end,4:(end-1)) = num2cell(Volume);
output(2:end,end) = num2cell(1+ceil(sum(Volume,2)*10)/10);
            
    
%%
cell2tsv(filename, output);


end

