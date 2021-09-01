 function varargout = ba(varargin)
%BA Bland-Altman Analysis of agreement between measurement methods.
%   Bland-Altman Analysis is a statistical method published by J. Martin
%   Bland and Douglas G. Altman in the 1980s and further developed later
%   on. It is used to determine the agreement between two measurement
%   methods that measure the same quantity. It is a popular method in
%   biostatistics and chemistry.
%
%   Syntax
%   s = BA(x,y) performs Bland-Altman Analysis on x and y, which are data
%   from two measurement methods for the same quantity respectively. x and
%   y must be numeric vectors of the same length. The calculations are done
%   at a significance level of alpha = 0.05. Output s is a structure
%   containing multiple fields with descriptive statistics about the
%   agreement of x and y. For more details on s see section Output below.
%
%   s = BA(x,y,alpha) specifies the significance level to calculate the
%   limits of agreement and confidence intervals with. alpha must be a
%   scalar in the interval [0,1]. If alpha is not specified a value of 0.05
%   is used by default to calculate 95% limits of agreement and confidence
%   intervals.
%
%   s = BA(__,Name,Value) specifies additional options using one or more
%   name-value pair arguments, in addition to any of the input arguments in
%   the previous syntaxes. For example, you can specify to create the
%   mean-difference plot using the 'PlotMeanDifference' name-value pair
%   argument.
%
%   BA(__) can be used to plot the data without returning an output
%   argument.
%
%   __ = BA(f,__) specifies the figure(s) f in which to create the plots
%   specified with the corresponding name-value pairs. The number of
%   figures in f must equal one or the number of specified plots.
%
%   __ = BA(ax,__) specifies the (array of) axes in which to create the
%   plots specified with the corresponding name-value pairs. The number of
%   axes in ax must equal the number of specified plots.
%
%   Examples
%   See and run the ba1999demo.m script for examples of the syntax of BA
%   used with data from the 1999 article by Bland and Altman. Calling BA
%   without input arguments also runs the demonstration script.
%
%   Name-Value Pair Arguments
%   Specify optional comma-separated pairs of Name,Value arguments to
%   access various options. Name is the argument name and Value is the
%   corresponding value. Name must appear inside single quotes (' '). You
%   can specify several name and value pair arguments in any order as
%   Name1,Value1,...,NameN,ValueN.
%   Example: 'XName','X','YName','Y'
%
%   'XName': Name of x variable
%   inputname of input argument x (default) | string
%   Name of x variable, specified as a string. 'XName' is used in the plot
%   titles.
%   Example: 'XName','X' sets the first measurement's name to 'X'.
%
%   'YName': Name of y variable
%   inputname of input argument y (default) | string
%   Name of y variable, specified as a string. 'YName' is used in the plot
%   titles.
%   Example: 'YName','Y' sets the second measurement's name to 'Y'.
%
%   'Exclude': Observation pairs to exclude
%   [] (default) | logical indices | numeric indices
%   Observation pairs to exclude, specified as logical or numeric indices
%   to index into x and y. The specified elements are removed from x and y
%   before any calculations or plots.
%   Example: 'Exclude',[1, 3, 4] excludes elements 1, 3 and 4 from x and y.
%   Example: 'Exclude',[0 0 1 0 1 1 0 0 1] excludes the true elements from
%   x and y.
%
%   'PlotMeanDifference': Create mean-difference plot
%   false (default) | true
%   Create the mean-difference plot if the specified value is true. The
%   mean-difference plot is a scatter plot of the difference between
%   observations versus their mean. Specifying the 'PlotAll' Name-Value
%   pair argument as true creates the mean-difference plot, regardless of
%   the 'PlotMeanDifference' value.
%
%   'PlotCorrelation': Create correlation plot
%   false (default) | true
%   Create the correlation plot if the specified value is true. The
%   correlation plot is a scatter plot of x and y. Specifying the 'PlotAll'
%   Name-Value pair argument as true creates the correlation plot,
%   regardless of the 'PlotCorrelation' value.
%
%   'PlotAll': Create all plots
%   false (default) | true
%   Create mean-difference and correlation plots if the specified value is
%   true. Setting 'PlotAll' to true overrides any value given to the
%   'PlotMeanDifference' and 'PlotCorrelation' Name-Value pair arguments.
%   However, setting it to false does not override the individual plot
%   Name-Value pair arguments.
%
%   'PlotStatistics': Add statistics to the created plots
%   'none' (default) | 'basic' | 'extended'
%   Add statistics to the created plots, specified as 'basic' or
%   'extended'. 'basic' specifies a basic set of statistics to add.
%   'extended' adds a more extended set of statistics to the plots. The
%   following statistics are added to the plots.
%   The basic set adds labelled lines for the limits of agreement to the
%   mean-difference plot. It also adds the Spearman rank correlation
%   coefficient to the legend of the mean difference plot. Furthermore, the
%   Pearson correlation coefficient is added to the legend of the
%   correlation plot.
%   The extended set adds the statistics of the basic set. Additionally,
%   the confidence intervals of the mean difference and limits of agreement
%   are plotted as error bars in the mean-difference plot. The extended set
%   does not add statistics other than the basic set to the correlation
%   plot.
%   If no plots are created, the 'PlotStatistics' value is ignored.
%
%   'PlotLeastSquares': Add a least-squares line to the created plots
%   false (default) | true
%   Add a least-squares line to the created plots if the specified value is
%   true. The least-sqares line describes the simple linear relationship
%   between x and y in the correlation plot and the mean and difference in
%   the mean-difference plot.
%
%   Output
%   The only output argument s is optional. It is a scalar structure
%   containing multiple fields with descriptive statistics about the
%   agreement of x and y. s contains the following fields:
%   
%   muD: the mean difference between x and y, also called the bias.
% 
%   muDCI: the 95% (default, depending on alpha) confidence interval of
%   the mean difference.
% 
%   loa: the 95% (default, depending on alpha) limits of agreement, a 2
%   element vector. The first element is the lower limit of agreement, the
%   second is the upper.
% 
%   loaCI: the 95% (default, depending on alpha) confidence interval of the
%   limits of agreement, a 2x2 matrix. The first column corresponds to
%   lower limit of agreement, the second to the upper limit. The first and
%   second row correspond to the lower and upper confidence interval bound
%   respectively.
% 
%   sD: the standard deviation of the differences.
% 
%   rSMuD: the Spearman rank correlation between mean and difference.
% 
%   pRSMuD: the p-value of the Spearman rank correlation for testing the
%   hypothesis of no correlation against the alternative that there is a
%   nonzero correlation.
%
%   About
%   This MATLAB function is an implementation of the methods in the 1999
%   article by Bland and Altman:
%   Bland & Altman 1999 - Measuring agreement methods in comparison studies
%   http://smm.sagepub.com/content/8/2/135.abstract
%
%   You might not have access to this article. Access it through your
%   institution's library or buy it.
%
%   The article comprises of 5 methodological sections (sections 2-6). The
%   current version of this MATLAB file (8 september 2016) implements the
%   first methodological section. More sections will be added in the
%   future.
%
%   The demonstration script `ba1999demo.m` is an implementation of the
%   calculations done by Bland and Altman in the article. Their article
%   contains a number of example data sets, which they use in their
%   methods. The demonstration script illustrates the same results and the
%   syntax used to obtain them.
%
%   See also BA1999DEMO

%   Copyright (C) 2016 Erik Huizinga, huizinga.erik@gmail.com
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% inputs
% demo if no input arguments
if ~nargin, ba1999demo, return, end

in = varargin;

% check if first input is (an array of) handles
if all(isgraphics(in{1}))
    h = in{1};
    in = in(2:end);
    ixname = 2;
    iyname = 3;
else
    h = [];
    ixname = 1;
    iyname = 2;
end

% prepare parser
p = inputParser;
p.addRequired('x',@isnumeric)
p.addRequired('y',@isnumeric)
p.addOptional('a',.05,@isnumeric)
p.addParameter('XName',inputname(ixname),@ischar)
p.addParameter('YName',inputname(iyname),@ischar)
p.addParameter('PlotAll',false,@validatelogical)
p.addParameter('PlotMeanDifference',false,@validatelogical)
p.addParameter('PlotCorrelation',false,@validatelogical)
p.addParameter('Exclude',[],@validatelogical)
p.addParameter('PlotStatistics','none',@ischar)
p.addParameter('PlotLeastSquares',false,@validatelogical)

% parse inputs
parse(p,in{:})
s2v(p.Results); %#ok<*NODEF>

%% validate and preprocess inputs
% x and y: measurements of two methods
if numel(x)~=numel(y)
    error 'Number of elements in x and y must be equal.'
end
% force column vectors
if isvector(x), x = x(:); end
if isvector(y), y = y(:); end

% alpha: significance level
if ~isscalar(a) && a<0 && a>1
    error 'alpha must be a scalar in the interval [0,1].'
end

% xName and yName
if ~iscellstr(XName), XName = cellstr(XName); end
if ~iscellstr(YName), YName = cellstr(YName); end
xName = strjoin(XName,', ');
yName = strjoin(YName,', ');

% validate plot arguments
[doPlotMD,axMD,doPlotC,axC] = validatePlotArgs( ...
    PlotAll, PlotMeanDifference, PlotCorrelation, h ...
    );

% exclude samples
lex = false(size(x));
lex(Exclude) = true;
x(lex) = [];
y(lex) = [];

% statistics set to plot
switch lower(PlotStatistics)
    case 'none'
        doPlotBasicStats = false;
        doPlotExtStats = false;
    case 'basic'
        doPlotBasicStats = true;
        doPlotExtStats = false;
    case 'extended'
        doPlotBasicStats = true;
        doPlotExtStats = true;
    otherwise
        error 'Unknown value for the ''PlotStatistics'' name-value pair.'
end

% least-squares line
doPlotLS = logical(PlotLeastSquares);

%% Bland-Altman analysis
out = baloa( ...
    x, xName, y, yName, a, ...
    doPlotMD, axMD, ...
    doPlotC, axC, ...
    doPlotBasicStats, doPlotExtStats, doPlotLS);

%% output
if nargout, varargout = out; end
end