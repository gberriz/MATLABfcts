function parts = regexpcellsplit(strcell, express, selected)
% pos = regexpcellsplit(strcell, expression, selected)
%   find the parts in a cell of strings strcell split based on the regular
%   expression express. Return the cell array parts of the same length as
%   strcell containing the parts of expression in strcell{i} or the whole string
%   if the expression doesn't match strcell. If multiple
%   parts are found, all are passed unless selected is set to a integer
%   (e.g. 1 of only the first part of the split is reported).
%
%   note: if selected is negative, works as -selected, BUT
%       and if -selected > length( splitted str), parts will be empty
%

parts = cell(size(strcell));
temp = regexp(strcell(:,1),express,'split');
for i=1:length(temp)
    if exist('selected','var')
        if selected<0 && -selected>length(temp{i})
            parts{i} = '';
        else
            parts{i} = temp{i}{min(abs(selected), length(temp{i}))};
        end
    elseif length(temp{i})==1
        parts{i} = temp{i}{1};
    else
        parts{i} = temp{i};
    end
end