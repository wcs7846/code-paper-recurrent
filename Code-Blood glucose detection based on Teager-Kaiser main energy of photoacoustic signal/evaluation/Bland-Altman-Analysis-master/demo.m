%BA1999DEMO Demonstration of:
% Bland and Altman 1999
% Measuring agreement methods in comparison studies
% http://smm.sagepub.com/content/8/2/135.abstract

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

%% display demo text
disp 'Demonstration of Bland-Altman analysis'
close all; clear; clc;
%% load data
addpath('./lib');
% Table 1 from the article
load fitting_data

JName = 'Real concentration of blood glucose';
SName = 'Predicted concentration of blood glucose';

%% 2 Limits of agreement (article page 136)
% This section only considers the first observations by J and S.
J1 = concentration_list';
S1 = pred_vec;

s = ba(J1,S1);
muD = s.muD; % mean difference (bias)
sD = s.sD; % standard deviation of difference
loa = s.loa; % limits of agreement
display(muD) % article: -16.29 mmHg
display(sD) % article: 19.61 mmHg
display(loa) % article: [-54.7, 22.1] mmHg

% exclude subjects 78 and 80 (article page number 139)
fprintf '\nExcluding subjects 78 and 80 (outliers):\n'
% iEx = [78,80];
% lEx = (n==78) | (n==80); % alternative: logical indices
% s = ba(J1,S1, 'Exclude',iEx);
% muD = s.muD;
% loa = s.loa;
% display(muD) % article: -14.9 mmHg
% Erratum: a small error, probably a typo, exist in the original article.
% Note the difference between the muD and the article's mean difference.
% However, because the limits of agreement are the same and, of course, the
% mean difference equals the mean of the limits of agreement, thus the
% calculation here appears to be correct.
display(loa) % article: [-43.6, 15.0] mmHg

%% 2.1 Graphical presentation of agreement (article page 140)
% If any other figures exist, the numbers below might not be correct. This
% script was written with all other figures closed.
% Figure 1 corresponds to article figure 1, figure 2 corresponds to article
% figure 2.
s = ba(figure, J1,S1, 'XName',JName, 'YName',SName, 'PlotAll',true);
rSMuD = s.rSMuD;
% Article: 0.07 (p. 140), here -0.03, so an error exists either here or in
% the article.
display(rSMuD); % Spearman rank correlation between mean and difference

% Figure 3 is equal to figure 2, but has added limits of agreement. It
% corresponds to article figure 3.
ba(figure, J1,S1, 'XName',JName, 'YName',SName, ...
    'PlotMeanDifference',true, 'PlotStatistics','basic')

%% 2.2 Precision of the estimated limits of agreement (article page 141)
clf % Clear figure 3 and recreate it with additional confidence intervals
% (not shown in article).
s = ba(gcf, J1,S1, 'XName',JName, 'YName',SName, ...
    'PlotMeanDifference',true, 'PlotStatistics','extended');
muDCI = s.muDCI;
loaCI = s.loaCI; % confidence interval of the loa
% article loaCI (p. 142):
% [-61.9, 14.9
%  -47.5, 29.3]
display(loaCI)
% article muDCI (p. 142): [-20.5 -12.1]
display(muDCI)