

fprintf('\n\nExample of importing the data:\n---------------------------\n\n');

fprintf('folder = ''./results/''; \n')

fprintf('t_data = Import_PlatesCellCountData([folder ''ColumbusOutput.txt''], ... \n');
fprintf('\t[folder ''PlateInfo.tsv''], ''NobjField'', ... %% example with multiple non-default fields \n');
fprintf('\t{''Nuclei_Hoechst_NumberOfObjects'' ''Nuclei_LDR_NumberOfObjects'' ''Nuclei_Hoechst_LDRpos_NumberOfObjects''});\n');

fprintf('\n%%%%\nt_annotated = Annotate_CellCountData(t_data, folder); %% where the treatment design files are stored \n');
fprintf('t_corrected = EdgeCorrecting_CellCountData(t_annotated); %% if edge correction is needed \n');

%
fprintf('\n%%%%\n[t_mean, t_processed] = Merge_CellCountData(t_corrected); \n');