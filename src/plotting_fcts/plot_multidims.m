function a = plot_multidims(t_data, varargin)

%  a = plot_multidims(t_data, varargin)
%       inputs: xplotkey, yplotkey, xaxiskey, yaxiskey, colorkey,
%       xtransform, ytransform, axischanges, mean_SEM, xspacing, yspacing,
%       yval_lines
%

Generate_Plotting_parameters

noyplot = 'noyplot';
nocplot = 'nocplot';

p = inputParser;
addParameter(p, 'xplotkey', @isstr)
addParameter(p, 'yplotkey', noyplot, @isstr)
addParameter(p, 'xaxiskey', @isstr)
addParameter(p, 'yaxiskey', @isstr)
addParameter(p, 'colorkey', nocplot, @isstr)
addParameter(p, 'xtransform', @(x)x, @(x)isa(x,'function_handle'));
addParameter(p, 'ytransform', @(x)x, @(x)isa(x,'function_handle'))
addParameter(p, 'axischanges', @(x) set(x,'fontsize',6), @(x) isa(x,'function_handle'))
addParameter(p, 'mean_SEM', true, @islogical)
addParameter(p, 'plotcolors', Plotting_parameters.colors, @isnumeric)
addParameter(p, 'xspacing', .03, @isnumeric)
addParameter(p, 'yspacing', .07, @isnumeric)
addParameter(p, 'yval_lines', [0 1], @isnumeric)

parse(p,varargin{:})
p = p.Results;

xplotkeys = unique(t_data.(p.xplotkey));
if strcmp(p.yplotkey, noyplot)
    yplotkeys = 1;
    t_data = [t_data table(ones(height(t_data),1), 'variablenames', {noyplot})];
else
    yplotkeys = unique(t_data.(p.yplotkey));
end
if strcmp(p.colorkey, nocplot)
    colorkeys = 1;
    t_data = [t_data table(ones(height(t_data),1), 'variablenames', {nocplot})];
else
    colorkeys = unique(t_data.(p.colorkey));
end

%
t_data = TableToCategorical(t_data, {p.xplotkey, p.yplotkey, p.colorkey});

%%

nCols = length(xplotkeys);
nRows = length(yplotkeys);

xspacing = p.xspacing;
axis_width = (1-(nCols+1)*xspacing)/nCols;
yspacing = p.yspacing;
axis_height = (1-(nRows+1)*yspacing)/nRows;

for ixp=1:nCols
    for iyp = 1:nRows
        
        a(iyp, ixp) = get_newaxes([xspacing*1.5+(ixp-1)*(axis_width+xspacing) ...
            yspacing*1.5+(nRows-iyp)*(axis_height+yspacing) axis_width axis_height],1,...
            'fontsize',6);
        
        title([p.xplotkey '=' AnyToString(xplotkeys(ixp)) '; ' ...
            p.yplotkey '=' AnyToString(yplotkeys(iyp))])

        xvals = t_data.(p.xaxiskey)(t_data.(p.xplotkey)==xplotkeys(ixp) & ...
            t_data.(p.yplotkey)==yplotkeys(iyp));
        xvals = p.xtransform(xvals);
        for i=1:length(p.yval_lines)
            plot([min(xvals(:)) max(xvals(:))], p.yval_lines(i)*[1 1], ...
                '-', 'color', [.8 .8 .8])
        end
                
        for iC = 1:length(colorkeys)
            subt = t_data(t_data.(p.xplotkey)==xplotkeys(ixp) & ...
                t_data.(p.yplotkey)==yplotkeys(iyp) & ...
                t_data.(p.colorkey)==colorkeys(iC),{p.xaxiskey p.yaxiskey});
            if isempty(subt)
                continue
            end

            subt.(p.xaxiskey) = p.xtransform(subt.(p.xaxiskey));
            subt.(p.yaxiskey) = p.ytransform(subt.(p.yaxiskey));
            subt = sortrows(subt, p.xaxiskey);

            if ~p.mean_SEM
            h(iC) = plot(subt.(p.xaxiskey), subt.(p.yaxiskey), '.-', ...
                'color', p.plotcolors(mod(iC-1,size(p.plotcolors,1))+1,:));
            else
                y_mean = collapse(subt, @mean, 'keyvars', {p.xaxiskey}, 'valvars', {p.yaxiskey});
                y_SEM = collapse(subt, @SEM,  'keyvars', {p.xaxiskey}, 'valvars', {p.yaxiskey});
                if height(y_SEM) == height(subt)
                    h(iC) = plot(subt.(p.xaxiskey), subt.(p.yaxiskey), '.-', ...
                        'color', p.plotcolors(mod(iC-1,size(p.plotcolors,1))+1,:));
                else
                    h(iC) = errorbar(y_mean.(p.xaxiskey), y_mean.(p.yaxiskey),...
                        y_SEM.(p.yaxiskey),'.-', ...
                        'color', p.plotcolors(mod(iC-1,size(p.plotcolors,1))+1,:));
                end
            end
        end

        p.axischanges(gca)
        if iyp==nRows, xlabel(gca,p.xaxiskey,'fontweight','bold'), end
        if ixp==1, ylabel(gca, p.yaxiskey,'fontweight','bold'), end

        if ixp==nCols && iyp==nRows
            hl = legend(h, strcat(p.colorkey, '=', AnyToString(colorkeys)), ...
                'fontsize',6, 'orientation', 'horizontal');
            set(hl, 'position', [.01 .005 .98 .03])
        end
    end
end

if nargout==0
 clear a
end