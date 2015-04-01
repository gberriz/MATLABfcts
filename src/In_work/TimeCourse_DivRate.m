
t_data = t_annotated;


trace_vars = setdiff(varnames(t_data),...
    {'Time' 'Date' 'Cellcount' 'RelCellCnt' 'RelGrowth' 'Ctrlcount'});
t_traces = unique(t_data(:,trace_vars));

t_rate = table();
for it = 1:height(t_traces)
    
    t_temp = sortrows(t_data(eqtable(t_data(:,trace_vars), t_traces(it,:)),:),'Time');
    
    Time = [0; mean([t_temp.Time(1:(end-1)) t_temp.Time(2:end)],2)];
    Cellcount = [t_temp.Cellcount(1)
        mean([t_temp.Cellcount(1:(end-1)) t_temp.Cellcount(2:end)],2)];
    dx = diff(t_temp.Cellcount);
    dt = diff(t_temp.Time)*24;
    DivRate = smooth([1e-4;dx./dt],7);
    
    t_rate = [t_rate;
        t_temp(:, trace_vars) table(Time, Cellcount, DivRate)];
end

ctrl_vars = ['Time'  ...
    setdiff(intersect(varnames(t_rate), trace_vars), ...
    {'Well' 'DrugName' 'HMSLid' 'Conc' 'pert_type' 'Untrt' 'TreatmentFile'})];
t_ctrl = collapse(t_rate(t_rate.Conc==0, ['DivRate' ctrl_vars]), @(x)mean(max(x,0)), ...
    'keyvars', ctrl_vars);

t_rate.RelDivRate = NaN(height(t_rate),1);
for i = 1:height(t_ctrl)
    idx = eqtable(t_ctrl(i,ctrl_vars), t_rate(:,ctrl_vars));
    t_rate.RelDivRate( idx ) = t_rate.DivRate( idx )./t_ctrl.DivRate(i);
end
t_rate.RelDivRate = max(min(t_rate.RelDivRate, 4), -2);

%%

get_newfigure(1,[50 50 1200 900])
plot_multidims(t_rate(t_rate.DrugName~='-',:), ...
    'xplotkey', 'SeedingDensity', 'yplotkey','DrugName',...
    'xaxiskey','Time', 'yaxiskey','RelDivRate','colorkey','Conc',...
    'axischanges',@(x)ylim(x,[-1 3]));

get_newfigure(2,[50 50 1200 900])
plot_multidims(t_rate(t_rate.DrugName~='-' & t_rate.Barcode=='140516_144932-V',:), ...
    'xplotkey', 'SeedingDensity', 'yplotkey','DrugName',...
    'xaxiskey','Time', 'yaxiskey','RelDivRate','colorkey','Conc',...
    'axischanges',@(x)ylim(x,[-1 3]));

get_newfigure(4,[50 50 1200 900])
plot_multidims(t_rate(t_rate.Barcode=='140516_144932-V',:), ...
    'xplotkey', 'SeedingDensity', 'yplotkey','DrugName',...
    'xaxiskey','Time', 'yaxiskey','DivRate','colorkey','Conc',...
    'axischanges',@(x)ylim(x,[-1 3]));
