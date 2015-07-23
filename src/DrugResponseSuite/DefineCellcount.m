function t_raw = DefineCellcount(t_raw, NobjField, p)

if ~isfield(p, 'Cellcount') || isempty(p.Cellcount)
    t_raw.Properties.VariableNames{NobjField{1}} = 'Cellcount';
elseif ~strcmp(p.Cellcount, 'none')
    temp = table2array(t_raw(:,NobjField));
    temp = p.Cellcount(temp);
    t_raw.Cellcount = temp;
end
