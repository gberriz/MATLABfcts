function lsgit(idx)
% lsgit(idx)

HomeFolder = fileparts(mfilename('fullpath'));

folders = dir(HomeFolder);
folders = setdiff({folders([folders.isdir]).name}', {'.' '..'});

fprintf('Folder %s:\n', HomeFolder);
fprintf([repmat('-',1,8+length(HomeFolder)) '\n']);

if ~exist('idx','var')
    for i=1:length(folders)
        fprintf('\n%2i: %s\n',i, folders{i});
        filenames = dir([HomeFolder filesep folders{i} '/*.m']);
        filenames = {filenames.name}';
        for j=1:length(filenames)
            fprintf('\t%-5.2f: %s\n',i+j/100, filenames{j});
        end
    end
elseif idx==0
    for i=1:length(folders)
        fprintf('\n%2i: %s\n',i, folders{i});
    end
else
    filenames = dir([HomeFolder filesep folders{floor(idx)} '/*.m']);
    
    fprintf('\nFolder %s (%i)\n', folders{floor(idx)}, floor(idx));
    filenames = {filenames.name}';
    if floor(idx)==idx
        for j=1:length(filenames)
            fprintf('* %s (%.2f)\n', filenames{j}, idx+j/100);
            a = help(filenames{j});
            fprintf(['\t' a(1:find(a==sprintf('\n'),1,'first')-1) '\n\n']);
        end
    else
        j = round(100*mod(idx,1));
        fprintf('* %s (%.2f)\n', filenames{j}, idx);
        help(filenames{j})
    end
end

