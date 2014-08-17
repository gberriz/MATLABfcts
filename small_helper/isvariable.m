function isVar = isvariable(t, VariableName)
isVar = ismember(VariableName, t.Properties.VariableNames);