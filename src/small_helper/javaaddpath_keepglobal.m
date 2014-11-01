function javaaddpath_keepglobal(dpath, varargin)
% javaaddpath_keepglobal(dpath, varargin)
%  same function as JAVAADDPATH DIRJAR but keep the global variables. It
%  adds the specified directory or jar file to the current dynamic Java
%  path.
%
%  When loading Java classes, MATLAB always searches the static Java path
%  before the dynamic Java path.  The static path is fixed at startup and
%  cannot be changed.  It contains, in the following order:
% 
%     - MATLAB's built-in Java path
%     - the contents of javaclasspath.txt in the startup directory
%     - the contents of javaclasspath.txt in the preferences directory
%       (see 'prefdir')
% 
%  Enter 'javaclasspath' to see the static and current dynamic paths.
%
%  JAVAADDPATH(DIRJAR)  ... adds directories or jar files
%  to the beginning of the current dynamic Java path.
%  Relative paths are converted to absolute paths.
%
%  JAVAADDPATH(..., '-END') appends the specified directories.
%
%  Use the functional form of JAVAADDPATH, such as 
%  JAVAADDPATH({'dirjar','dirjar',...}), when the directory 
%  specification is stored in a string or cell array of
%  strings.
%

g = who('global');
BaseVars = evalin('base','who');
for i=1:length(g)
    eval(sprintf('global %s', g{i}))
    eval(sprintf('%s_x = %s;', g{i}, g{i}))
end

javaaddpath(dpath, varargin{:})

for i=1:length(g)
    eval(sprintf('global %s', g{i}))
    eval(sprintf('%s = %s_x;', g{i}, g{i}))
    if ismember(g{i}, BaseVars)
        eval(sprintf('assignin(''base'', g{i}, %s);', g{i}))
    end
end
