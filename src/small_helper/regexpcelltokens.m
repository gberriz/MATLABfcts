function toks = regexpcelltokens(strcell, express, selected)
% pos = regexpcelltokens(strcell, expression, selected)
%   find the token defined by the regular expression express in a cell
%   of strings strcell. Return the cell array toks of the same length as
%   strcell containing the token (if any) of expression in strcell{i} or ''
%   (empty string) if the expression doesn't match strcell. If multiple
%   tokens are found, all are passed unless selected is set to a integer
%   (e.g. 1 of only the first token match needs to be reported).
%
%

toks = cell(size(strcell));
temp = regexp(strcell(:,1),express,'tokens');
for i=1:length(temp)
    if isempty(temp{i}) || all(cellfun(@isempty,temp{i}))
        toks{i} = '';
    elseif length(temp{i})==1
        toks(i) = temp{i}{1};
    elseif exist('selected','var')
        toks(i) = temp{i}{min(selected, length(temp{i}))};
    else
        toks{i} = temp{i};
    end
end