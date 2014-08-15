function cstr = table2cellstr(t_in, headers)

% for i=1:size(a,2)
%     if iscategorical(a.(a.Properties.VariableNames{i})(1))
%         a.(a.Properties.VariableNames{i}) = cellstr(a.(a.Properties.VariableNames{i}));
%     elseif isnumeric(a.(a.Properties.VariableNames{i})(1))
%         a.(a.Properties.VariableNames{i}) = cellfun_brdcast(@num2str,...
%             num2cell(a.(a.Properties.VariableNames{i})));
%     end
% end
t_in = TableToString(t_in, 0);

if ~exist('headers','var') || headers
    cstr = [t_in.Properties.VariableNames;    table2cell(t_in)];
else
    cstr = table2cell(t_in);
end
