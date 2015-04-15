function t_rate = TimeCourse_DivRate(t_data, plate_inkeys)
% t_rate = TimeCourse_DivRate(t_data, plate_inkeys)
%
%   process the data from cell count. The data will be split for controls according
%   to plate_keys (with 'Barcode' 'CellLine' 'Time' as mandatory and
%   default).
%   Needs also the column 'pert_type' (looking for value
%   'ctl_vehicle' or 'trt_cp'/'trt_poscon'), 'Barcode' and 'Well' as unique
%   identifiers
%
%   outputs are:    t_rate with the values for all treated wells
%

% unique identifiers for following timecourse
trace_vars = {'Barcode' 'Well'};

if exist('plate_inkeys','var') && ~isempty(plate_inkeys)
    plate_keys = unique([{'CellLine' 'Barcode' 'Time'} plate_inkeys]);
else
    plate_keys = {'CellLine' 'Barcode' 'Time'};
    plate_inkeys = {};
end
%%

t_location = unique(t_data(:,trace_vars));

t_rate = table();
for it = 1:height(t_location)
    % get all time points for each well
    t_temp = sortrows(t_data(eqtable(t_data(:,trace_vars), t_location(it,:)),:),'Time');
    
    Time = [0; mean([t_temp.Time(1:(end-1)) t_temp.Time(2:end)],2)];
    Cellcount = [t_temp.Cellcount(1)
        mean([t_temp.Cellcount(1:(end-1)) t_temp.Cellcount(2:end)],2)];
    dx = diff(t_temp.Cellcount);
    dt = diff(t_temp.Time)/24;
    DivRate = dx./dt./Cellcount(2:end);
    DivRate = smooth([mean(DivRate(1:3));DivRate], 3);

    n = NaN(width(t_temp),1);
    for i=1:width(t_temp), n(i) = length(unique(t_temp.(i))); end;
    annotation_vars = setdiff(t_temp.Properties.VariableNames(n==1), ...
        [trace_vars 'RelCellCnt' 'RelGrowth']);

    t_rate = [t_rate;
        t_temp(:, trace_vars) table(Time, Cellcount, DivRate) t_temp(:,annotation_vars)];
end

%
t_ctrl = collapse(t_rate(t_rate.pert_type=='ctl_vehicle', ['DivRate' plate_keys]), ...
    @(x)mean(max(x,0)), 'keyvars', plate_keys);

t_rate.RelDivRate = NaN(height(t_rate),1);
for i = 1:height(t_ctrl)
    idx = eqtable(t_ctrl(i,plate_keys), t_rate(:,plate_keys));
    t_rate.RelDivRate( idx ) = t_rate.DivRate( idx )./t_ctrl.DivRate(i);
end
t_rate.RelDivRate = max(min(t_rate.RelDivRate, 4), -2);


