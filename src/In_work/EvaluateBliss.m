function [BlissScore, Bliss, Results] = EvaluateBliss(t_data, varname, BlissType)
% [BlissScore, Bliss] = EvaluateBliss(t_data, varname, BlissType)
%
% Inputs:
%   t_data:     table with the columns: DrugName, Conc, DrugName2, Conc2,
%                   RelGrowth/(varname)
%               (varname) should be either increasing with dose or
%                   normalized to 1 (control)
%
%   varname:    Variable on which to apply the Bliss evaluation (default =
%                   RelGrowth)
%
%   BlissType:  1 if normalized to one --> Bliss as (1-x)+(1-y) (default)
%               0 if not normalized --> Bliss as x+y
%
% Outputs:
%   BlissScore: average Bliss excess across all doses
%   Bliss:      matrix of Bliss excess
%   Results:    matrix for the values of RelGrowth/(varname)
%
% Improvment to make:
%   - support other input formats
%   - currently assuming that all doses used in the combo are also in the
%       single agent
%
%

if ~exist('varname','var') || isempty(varname)
    varname = 'RelGrowth';
else
    assert(isvariable(t_data,varname), ...
        '%s is not a header of the table t_data', varname)
end

if isvariable(t_data,'CellLine')
    assert(length(unique(t_data.CellLine))==1, 'There are more than one cell line: %s', ...
        strjoin(cellstr(unique(t_data.CellLine)), ', '))
end

Drugs = unique(t_data.DrugName(t_data.DrugName2~='-'));
assert(~isempty(Drugs), 'No ''Drugs(1)'' found')
assert(length(Drugs)==1, 'There is more than one ''Drug1'': %s', ...
    strjoin(cellstr(Drugs),', '))

temp = setdiff(unique(t_data.DrugName(t_data.DrugName2=='-')), Drugs);
assert(length(temp)==1, 'There is more than one ''Drug2'': %s', ...
    strjoin(cellstr(temp),', '))
Drugs(2) = temp;

%%
Conc1 = unique([t_data.Conc(t_data.DrugName==Drugs(1))
    t_data.Conc2(t_data.DrugName2==Drugs(1))]);
Conc2 = unique([t_data.Conc(t_data.DrugName==Drugs(2))
    t_data.Conc2(t_data.DrugName2==Drugs(2))]);
Concs = [{Conc1} {Conc2}];


Results = NaN(length(Conc1)+1, length(Conc2)+1);

temp = t_data(t_data.DrugName==Drugs(1) & t_data.DrugName2=='-',:);
temp = sortrows(temp , 'Conc');
[~,Cidx] = ismember(temp.Conc, Conc1);
Results(1+Cidx,1) = temp.(varname);


temp = t_data(t_data.DrugName==Drugs(2) & t_data.DrugName2=='-',:);
temp = sortrows(temp , 'Conc');
[~,Cidx] = ismember(temp.Conc, Conc2);
Results(1,1+Cidx) = temp.(varname);


t_sub = t_data(t_data.DrugName==Drugs(1) & t_data.DrugName2==Drugs(2),:);
Concs = unique(t_sub.Conc2)';
for iDo = 1:length(Concs)
    temp = sortrows(t_sub(t_sub.Conc2 == Concs(iDo),:),'Conc');
    [~,Cidx] = ismember(temp.Conc, Conc1);
    Results(1+Cidx, 1+find(Conc2==Concs(iDo))) = temp.(varname);
end

%%
if (mean(diff(Results(:,1))>0)>.5 && mean(diff(Results(:,1))>0)>.5 && ...
        mean(Results(:)>0)>.5 && ~strcmp(varname,'RelGrowth')) || ...
        (exist('BlissType','var') && BlissType==0)
    % generally increasing with dose --> calculate Bliss as x*y
    Bliss = Results - repmat(Results(1,:),size(Results,1),1) - ...
        repmat(Results(:,1),1,size(Results,2)) + ...
        repmat(Results(1,:),size(Results,1),1).*repmat(Results(:,1),1,size(Results,2));
    Results(1,1) = 0;
elseif mean(Results(:)<1)>.5 || any(ismember(varname,{'RelGrowth' 'RelCellCnt'})) || ...
        (exist('BlissType','var') && BlissType==1)
    % mainly decreasing with dose, normalized at 1 --> calculate Bliss as (1-x)*(1-y)
    % of Relative cell growth/count
    Bliss = 1-Results - (1-repmat(Results(1,:),size(Results,1),1)) - ...
        (1-repmat(Results(:,1),1,size(Results,2))) + ...
        (1-repmat(Results(1,:),size(Results,1),1)).*(1-repmat(Results(:,1),1,size(Results,2)));
    Results(1,1) = 1;
else
    error('(varname) is either not increasing or not normalized')
end

temp = Bliss(2:end,2:end);
BlissScore = nanmean(temp(:));

