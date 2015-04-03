function merge_pdf(inputfiles, outputfile, args)
%  merge_pdf(inputfiles, outputfile, args)
%
%   relies on the 'server' version of pdftk (see
%   https://www.pdflabs.com/tools/pdftk-server/) pdftk should be accesible
%   from the current folder
%
%   inputfiles is a string readable in DOS formet (like ' any*.pdf ')
%   outpfile is the name of a pfd file for the output
%   args is optional arguments for the pdftk command (default is ' cat ')
%
%   WARNING: overwritting the output file without notice
%

if ~exist('args','var')
    args = 'cat';
end




d = pwd;
% absolute path needed for MS-DOS
inputfiles = [d filesep ReplaceName(inputfiles, '/\', filesep)];
outputfile = [d filesep ReplaceName(outputfile, '/\', filesep)];

cmd = sprintf('pdftk %s %s output %s', inputfiles, args, outputfile);

status = system(cmd);
assert(status==0, 'pdftk was not successful')
