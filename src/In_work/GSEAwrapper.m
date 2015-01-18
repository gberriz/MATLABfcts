
function [t_results, Genelist] = GSEAwrapper(Genelist, Geneset, varargin)
% [t_results, Genelist] = GSEAwrapper(Genelist, Geneset, varargin)
%
%   Inputs: 
%   - Genelist:
%       - filename
%       - cell array or table (first column is gene symbol,
%               2nd column is weight)
%   - Geneset:
%       - tag (GObp, GOcc, GOmf, KEGG, Reactome)
%       - filename
%       - cell array or table with genes
%   - options:
%       - Outputfolder (folder in which the GSEA output folder will be
%           saved default is none, results are not save) 
%       - Nplot (number of sets to display for html output) 
%       - label (for the GSEA output, default is geneset label)
%       - Randseed (value for randomization)
%
%   Output:
%   - t_results: table with results for the selected sets
%   - Genelist : ordered gene lists
%
% requires the following files from GSEA to be stored in a unique folder:
%   - gsea2-2.1.0.jar
%   - gene set list such as (default), that can be called with [xx]
%       - c5.mf.v4.0.symbols.gmt [GOmf]
%       - c5.cc.v4.0.symbols.gmt [GOcc]
%       - c5.bf.v4.0.symbols.gmt [GObf]
%       - c2.cp.kegg.v4.0.symbols.gmt [KEGG]
%       - c2.cp.reactome.v4.0.symbols.gmt [Reactome]
%
% For more options, see: 
% http://www.broadinstitute.org/gsea/doc/GSEAUserGuideTEXT.htm#_Running_GSEA_from

global WorkFolder

% default values (can be specified as input parameters)
p = inputParser;
addParameter(p,'Randseed',ceil(mod(now,1)*1e5), @isscalar);
addParameter(p,'nperm',1e3, @isscalar);
addParameter(p,'Nplot',0, @isscalar);
addParameter(p,'set_min',15, @isscalar);
addParameter(p,'set_max',500, @isscalar);
addParameter(p,'Outputfolder', '', @ischar);
addParameter(p,'label', '', @ischar);
addParameter(p,'verbatim', false, @islogical);
addParameter(p,'GSEAfolder', [WorkFolder 'GSEA_java' filesep], @ischar);

parse(p,varargin{:});
p = p.Results;
Outputfolder = p.Outputfolder;
label = p.label;
GSEAfolder = p.GSEAfolder;

% create a temp folder based on the day (will be anyway created by GSEA
tempfolder = ['.' filesep ...
    lower(char(datetime('now','format','MMMd'))) filesep];
if ~exist(tempfolder,'dir')
    mkdir(tempfolder)
else
    delete([tempfolder filesep '*'])
end
if isempty(Outputfolder)
    Outputfolder = tempfolder;
    Nplot = 0;
elseif ~exist(Outputfolder,'dir')
    mkdir(Outputfolder)
end

Inputfile = [tempfolder 'temp.rnk'];
if ischar(Genelist)
    assert(exist(Genelist,'file')>0, 'Input file not found')
    Inputfile = Genelist;
elseif istable(Genelist)
    table2tsv(Genelist, Inputfile, 0);
elseif iscellstr(Input)
    cell2tsv(Inputfile, Genelist);
else
    error('unknown format for Input')
end


Setfile = [tempfolder 'tempset.gmt'];
switch Geneset
    case 'GOmf'
        Setfile = [GSEAfolder 'c5.mf.v4.0.symbols.gmt'];
    case 'GObp'
        Setfile = [GSEAfolder 'c5.bp.v4.0.symbols.gmt'];
    case 'GOcc'
        Setfile = [GSEAfolder 'c5.cc.v4.0.symbols.gmt'];
    case 'KEGG'
        Setfile = [GSEAfolder 'c2.cp.kegg.v4.0.symbols.gmt'];
    case 'Reactome'
        Setfile = [GSEAfolder 'c2.cp.reactome.v4.0.symbols.gmt'];
    otherwise
        if ischar(Geneset)
            assert(exist(Geneset,'file'), 'Geneset file not found')
            [~,~,ext] = filepart(Geneset);
            assert(strcmp(ext,'gmt'),'Geneset file must be a .gmt format')
            Setfile = Geneset;
        elseif istable(Geneset)
            table2tsv(Geneset, Setfile);
        elseif iscellstr(Geneset)
            cell2tsv(Setfile, Geneset);
        else
            error('unknown format for Geneset')
        end
end
% check if the file exist locally, otherwise use the
% Broad http address for the file: ftp.broadinstitute.org://pub/gsea/gene_sets/

% add the dataset in the label if not present
if ismember(Geneset, {'GOmf' 'GObp'})
    if isempty(strfind(label, Geneset))
        label = [label '_' Geneset];
    end
    if ~exist(Setfile, 'file')
        [~,Setfile] = fileparts(Setfile);
        Setfile = ['ftp.broadinstitute.org://pub/gsea/gene_sets/' Setfile '.gmt'];
    end
end
    
%% construct the command
cmd = [ ...
'java -Xmx512m -cp ' GSEAfolder 'gsea2-2.1.0.jar' ...
' xtools.gsea.GseaPreranked' ...
' -gmx ' Setfile ...
' -collapse false -mode Max_probe -norm meandiv -nperm ' num2str(p.nperm) ...
' -rnk ' Inputfile  ' -scoring_scheme weighted -rpt_label ' label ...
' -include_only_symbols true -make_sets true -plot_top_x ' num2str(p.Nplot) ...
' -rnd_seed ' num2str(p.Randseed) ' -set_max ' num2str(p.set_max) ...
' -set_min ' num2str(p.set_min) ' -zip_report false' ...
' -out ' Outputfolder ' -gui false'];

% and run it
if p.verbatim
    status = system(cmd);    
    disp('----------------------------------------------')
    assert(status==0, 'GSEA failed: %s')
else
    [status, output] = system(cmd);
    assert(status==0, 'GSEA failed: %s', output)
end
%% get the right folder

f = dir([Outputfolder filesep label '*']);
Lastrun = f(argmax(cellfun(@datenum,{f.date}))).name;
if isempty(Outputfolder)
    disp(['Results stored in ' Outputfolder])
    disp([' --->  ' filesep Lastrun])
else
    disp(' ---> GSEA done; cleaning stuff')
end

%% and get the genes of each set and their scores

% scores are in edb/results.edb
res = XMLElement.parse([Outputfolder filesep Lastrun filesep ...
    'edb' filesep 'results.edb']);

f = dir([Outputfolder filesep Lastrun filesep 'edb' filesep '*rnk']);
Genelist = tsv2cell([Outputfolder filesep Lastrun filesep 'edb' filesep f.name]);
%%
fields = {'GENESET' 'ES' 'NES' 'NP' 'FDR' 'FWER' 'HIT_INDICES' 'ES_PROFILE'};
fieldnames = {'GeneSet' 'Escore' 'NormEscore' 'pval' 'FDR' 'FWER' 'GeneIdx' ...
    'GeneEscore' 'GeneNames' 'LeadEdge'};

values = cell(length(res.children), length(fieldnames));
for i=1:length(res.children)
    for j=1:length(fields)
        values{i,j} = res.children(i).get(fields{j});
    end
    
    values{i,strcmp(fieldnames, 'GeneNames')} =  ...
        Genelist(values{i,strcmp(fields, {'HIT_INDICES'})},1)';
end

idx = ismember(fields, {'ES' 'NES' 'NP' 'FDR' 'FWER'});
values(:,idx) = cellfun2(@str2num,values(:,idx));

idx = strcmp(fields, 'GENESET');
values(:,idx) = regexpcellsplit(values(:,idx), '#', 2);
    
idx = strcmp(fields, 'ES_PROFILE');
values(:,idx) = cellfun2(@(x) cellfun(@str2num,regexp(x,' ','split')), ...
    values(:,idx));
idx = strcmp(fields, 'HIT_INDICES');
values(:,idx) = cellfun2(@(x) 1+cellfun(@str2num,regexp(x,' ','split')), ...
    values(:,idx));

t_results = cell2table(values, 'variablenames', fieldnames);

for i=1:height(t_results)
    t_results.GeneNames{i} = Genelist(t_results.GeneIdx{i},1)'; 
    if t_results.Escore(i)>0
        t_results.LeadEdge{i} = t_results.GeneNames{i}(1:find(...
            t_results.GeneEscore{i}==max(t_results.GeneEscore{i})));
    else
        t_results.LeadEdge{i} = t_results.GeneNames{i}(find(...
            t_results.GeneEscore{i}==min(t_results.GeneEscore{i})):end);
    end
end
% clean up
if exist(tempfolder,'dir')
    rmdir(tempfolder,'s')
end
        