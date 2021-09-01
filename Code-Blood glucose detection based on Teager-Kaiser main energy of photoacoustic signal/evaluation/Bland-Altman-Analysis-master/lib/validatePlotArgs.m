function [doPlotMD,axMD,doPlotC,axC] = validatePlotArgs( ...
    PlotAll, PlotMeanDifference, PlotCorrelation, h ...
    )
% default axes variables
axMD = [];
axC = [];

% doAllPlots
doPlotAll = logical(PlotAll);
if doPlotAll
    doPlotMD = true; % do mean-difference plot
    doPlotC = true; % do correlation plot
else
    doPlotMD = logical(PlotMeanDifference);
    doPlotC = logical(PlotCorrelation);
end

% validate number and type of handles in h for the requested plots
doPlot = [doPlotMD doPlotC];
if any(doPlot)
    nPlot = nnz(doPlot);
    if isempty(h)
        figure
        ax = gobjects(0);
    else
        nh = numel(h);
        if isaxes(h)
            if nh==nPlot
                ax = h;
            else
                error(['The number of requested plots (' num2str(nPlot) ...
                    ') is not equal to the number of axes (' ...
                    num2str(nh) ') in the first input argument.'])
            end
        else % h must be a figure or contain nplot figures
            if nh==1
                figure(h) % error if h is not a figure
                ax = gobjects(0);
            elseif nh==nPlot
                for n = nh:-1:1
                    figure(h(n))
                    ax(n) = axes;
                end
            elseif all(isfigure(h))
                error(['The number of requested plots (' num2str(nPlot) ...
                    ') is not equal to the number of figures (' ...
                    num2str(nh) ') in the first input argument.'])
            else
                error(['The first optional input argument must be an ' ...
                    '(array of) axes or figure handle(s).'])
            end
        end
    end
    
    % create axes if not yet existent
    if isempty(ax)
        if nPlot>1
            for n = nPlot:-1:1
                ax(n) = subplot(1,nPlot,n);
            end
        else
            ax = axes;
        end
    end
    
    % store axes in plot order
    if doPlotMD, axMD = ax(1); ax(1) = []; end
    if doPlotC, axC = ax(1); end % ax(1) = []; end
end
end