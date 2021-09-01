%% Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
% This script is designed to model the signal model of Blood glucose signal
% Author: MarkLHF
% Date: 2019-12-19(first version)
% Reference:

close all;
clear; format long;

% switch
debug_PropretyAnalysis = 1;
debug_RegressionAnalysis = 0;
debug_ConcentrationAnalysis = 1;

addpath('./lib');
addpath(genpath('./evaluation'));
% load the sensor data
load('real-data');
concentration_list = tmp_con'/10;

lineColor = linspecer(numel(concentration_list)); label=[];
% figure allocation
if debug_PropretyAnalysis
    show_rd = figure; % real data figure;
    show_feature1 = figure; % feature figure;
    show_feature2 = figure;
    show_feature3 = figure;
    show_feature4 = figure;
    n = 26;
end
tk_mat = [];
for n = length(concentration_list):-1:1
    tmpdata = data(n,:);
    sample_rate = 50e6; % 50MHz
    
    xrange = (1:numel(tmpdata))*1/sample_rate;
    if debug_PropretyAnalysis
        figure(show_rd);plot(xrange, tmpdata, 'color', lineColor(n, :));hold on;
        title(strcat('the real blood glucose signal: ',num2str(concentration_list(n)),'mg/dL'));
        xlabel('time/s');ylabel('pressure');
        
        label = [label, string(strcat(num2str(concentration_list(n)), ' mg/dL'))];
    end
    
    [tf_Stran, ts, fs] = st(tmpdata, 1);
    
    % frequency slice
    freq_vect = fs*sample_rate;
    
    upper_freq = 7e6; % 5MHz
    low_freq = 1.5e6; % 1MHz
    freq_band = find((freq_vect <= upper_freq) & (freq_vect >= low_freq));
    tf_Stran_band = tf_Stran(freq_band, :);
    
    % debug for time-frequence transfrom
    %     figure;imagesc(ts(freq_band)/sample_rate, fs(freq_band)*sample_rate, abs(tf_Stran_band));
    %     set(gca,'ydir','normal');colorbar;xlabel('Time');ylabel('Frequence');
    %     title(strcat('normal S transform ', num2str(concentration_list(n)), ' mg/dL'));
    
    tk = teagerEnergy(tf_Stran_band);
%     tk = teagerEnergy(tf_Stran_band);
    tk_mat = [tk_mat; tk];
    
    show_x = find(xrange > 000 & xrange < 5000 == 1) ;
    if debug_PropretyAnalysis
        figure(show_feature3); plot(xrange(show_x), tk(show_x),   'color', lineColor(n, :)); hold on;
    end
end

if debug_PropretyAnalysis
    figure(show_rd);
    colorbar('Ticks',0:1/(length(concentration_list)-1):1,'TickLabels', fliplr(label));colormap(lineColor);
    
    figure(show_feature3);title('Teager Energy');
    xlabel('time/s');ylabel('proprety value');xlim([min(xrange(show_x)),max(xrange(show_x))]);
    colorbar('Ticks',0:1/(length(concentration_list)-1):1,'TickLabels', fliplr(label));colormap(lineColor);
end
%% special location analysis
special_loc = 2.76e-6; % unit: um
loc = find(abs(xrange-special_loc) == min(abs(xrange-special_loc)));
if debug_ConcentrationAnalysis   
    figure;plot(fliplr(concentration_list),(tk_mat(:, loc)), 'ro');title('Teager Energy');hold on;
    xlabel('concentration(mg/dL)');ylabel('proprety value');%xlim([min(concentration_list),max(concentration_list)]);
   
end
%%
pp = loc;
% select tk as the data
data_select = (tk_mat(:, loc)); % data_select = x; concentration_list = y
data_select = data_select(length(data_select):-1:1);
% normalize - tk attributes
max_attr = max(data_select); min_attr = min(data_select);
% data_select = (data_select - min_attr)/(max_attr - min_attr);
% normalize - concentration
max_c = max(concentration_list); min_c = min(concentration_list);
% c = (concentration_list - min_c)/(max_c - min_c);
c = concentration_list;
    
% Mdl = fitlm(data_select, c);
% pred_vec = predict(Mdl, data_select);

data_ex = [data_select'; ones(1, length(data_select))]';
para_vec = inv(data_ex'*data_ex)*data_ex'*c';
Num_plot = 100;tt = min_attr:(max_attr-min_attr)/(Num_plot-1):max_attr;
line_vec = zeros(1, numel(tt));
for n = 1:numel(tt)
    tmp = [tt(n), 1]*para_vec;
    line_vec(n) = (tmp);
end
%plot(line_vec, tt, 'b-');hold on;
ylim([min(tt)*0.9, max(tt)*1.05]);
xlim([min(concentration_list)*0.5, max(concentration_list)*1.2]);
%legend('real data','fitting curve');

Endpoint_x = [min(concentration_list)*0.5, max(concentration_list)*1.2];
for n = 1:length(Endpoint_x)
    tmp = (Endpoint_x(n) - para_vec(2))/(para_vec(1) +eps);
    Endpoint_y(n) = tmp;
end
plot(Endpoint_x, Endpoint_y, 'b-');hold on;
legend('real data','fitting curve');
pred_vec = zeros(1, length(data_select));
for n = 1:length(data_select)
    tmp = [data_select(n), 1]*para_vec;
    pred_vec(n) = tmp;
end

% show
figure;
plot(1:length(concentration_list), concentration_list, 'r-');hold on;
plot(1:length(concentration_list), pred_vec, 'bo'); hold on;
legend('real', 'predict');
save('fitting_data.mat', 'pred_vec', 'concentration_list');

% 
mean_value = mean(pred_vec - concentration_list);
var_value  = std(pred_vec - concentration_list);
re = (pred_vec - concentration_list)./concentration_list;

mean_value_re = mean(re);
var_value_re = std(re);