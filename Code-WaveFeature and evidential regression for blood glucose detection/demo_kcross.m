%% Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
% This script is designed to model the signal model of Blood glucose signal
% Author: MarkLHF
% Date: 2019-12-19(first version)
% Reference:

close all;
clear; format long;

% switch
debug_PropretyAnalysis = 0;
debug_RegressionAnalysis = 1;
debug_ConcentrationAnalysis = 0;

addpath('./lib');
addpath(genpath('./evaluation'));

switch_data = 'long'; % chen or Long
if strcmpi(switch_data, 'chen')% load the sensor data - chen
    load('real-data');
    concentration_list = tmp_con'/10;
elseif strcmpi(switch_data, 'Long')% load the sensor data - Long
    load('real-data_total.mat'); % mg/dL
    concentration_list = tmp_con;
end

N = length(concentration_list);
lineColor = linspecer(numel(concentration_list)); label=[];
lineColor = lineColor(N:-1:1,:);
% figure allocation
if debug_PropretyAnalysis
    show_rd = figure; % real data figure;
    show_feature1 = figure;
%     show_feature2 = figure;
%     show_feature3 = figure;
end
feature_vec= [];
%% time-frequency features extraction
tk_mat = []; rmsf_mat = []; ibw_mat = []; fag_mat = []; loc_mat =[];
for n = 1:N % each sample
    tmpdata = data(n,:);
    sample_rate = 100e6; % 100MHz
    
    xrange = (1:numel(data(n,:)))*1/sample_rate;
    %     if debug_PropretyAnalysis
    %         figure(show_rd);plot(xrange, data(n,:), 'color', lineColor(n, :));hold on;
    %         title(strcat('the real blood glucose signal: ',num2str(concentration_list(n)),'mg/dL'));
    %         xlabel('time/s');ylabel('pressure');
    %
    %         label = [label, string(strcat(num2str(concentration_list(n)), ' mg/dL'))];
    %     end
    
    [tf_Stran, ts, fs] = st(data(n,:), 1);
    
    % frequency slice
    freq_vect = fs*sample_rate;
    
    upper_freq = 7e6; % 5MHz
    low_freq = 1.5e6; % 1MHz
    freq_band = find((freq_vect <= upper_freq) & (freq_vect >= low_freq));
    tf_Stran_band = tf_Stran(freq_band, :);
    
    % debug for time-frequence transfrom
    %         figure;imagesc(ts(freq_band)/sample_rate, fs(freq_band)*sample_rate, abs(tf_Stran_band));
    %         set(gca,'ydir','normal');colorbar;xlabel('Time');ylabel('Frequence');
    %         title(strcat('normal S transform ', num2str(concentration_list(n)), ' mg/dL'));
    
    tk = teagerEnergy(tf_Stran_band);
    fag = freqAttenGrad(tf_Stran_band, freq_vect(freq_band));
    [ cf, rmsf, ibw ] = centerFreq(tf_Stran_band, fs(freq_band)*sample_rate);
    
    tk_mat = [tk_mat; tk];
    rmsf_mat = [rmsf_mat; rmsf];
    ibw_mat = [ibw_mat; ibw];
    fag_mat = [fag_mat; fag];
    
    show_x = find(xrange > 000 & xrange < 5000 == 1) ;
    if debug_PropretyAnalysis
        figure(show_feature1); plot(xrange(show_x), tk(show_x),   'color', lineColor(n, :)); hold on; ylabel('TK energy');
%         figure(show_feature2); plot(xrange(show_x), rmsf_mat(show_x),   'color', lineColor(n, :)); hold on; ylabel('RMSF energy');
%         figure(show_feature3); plot(xrange(show_x), ibw_mat(show_x),   'color', lineColor(n, :)); hold on; ylabel('IBW energy');
    end
    loc_mat = [loc_mat; find(tk == max(tk))];
end
% time-frequency feature
% calculate the optimal location
[~,col] = size(tk_mat);

% use Pearson correlation coefficient
tk_coeff_vec = corr(concentration_list', tk_mat);
rmsf_coeff_vec = corr(concentration_list', rmsf_mat);
ibw_coeff_vec = corr(concentration_list', ibw_mat);

loc_tk = find(max(tk_coeff_vec) == tk_coeff_vec);
loc_rmsf = find(max(rmsf_coeff_vec) == rmsf_coeff_vec);
loc_ibw = find(max(ibw_coeff_vec) == ibw_coeff_vec);
% 使用时频特征作为特征向量
tf_vec = [tk_mat(:, loc_tk), rmsf_mat(:, loc_rmsf), ibw_mat(:, loc_ibw)]; %
feature_vec = [feature_vec tf_vec];
if debug_ConcentrationAnalysis
    figure;plot(concentration_list, tk_mat(:, loc), 'ro');title('Teager Energy');hold on;
    xlabel('concentration(mg/dL)');ylabel('proprety value');%xlim([min(concentration_list),max(concentration_list)]);
end
if debug_PropretyAnalysis
    figure(show_rd);
    colorbar('Ticks',0:1/(length(concentration_list)-1):1,'TickLabels', fliplr(label));colormap(lineColor);
    
    figure(show_feature3);title('Teager Energy');
    xlabel('time/s');ylabel('proprety value');xlim([min(xrange(show_x)),max(xrange(show_x))]);
    colorbar('Ticks',0:1/(length(concentration_list)-1):1,'TickLabels', fliplr(label));colormap(lineColor);
end
%% time domain features extraction
% Since there is no strong interference, the peak to peak value is equal to the difference between the maximum value and the minimum value.
MaxValue = max(data, [], 2);MinValue = min(data, [], 2);
ppk_loc_pmax = zeros(N,1); ppk_loc_pmin = zeros(N,1);
for n = 1:N % each sample
    [~,t_ppk_loc_pmax] = find(MaxValue(n) == data(n,:));
    [~,t_ppk_loc_pmin] = find(MinValue(n) == data(n,:));
    
    if t_ppk_loc_pmax > 1
        ppk_loc_pmax(n) = t_ppk_loc_pmax(1);
    end
    if t_ppk_loc_pmin > 1
        ppk_loc_pmin(n) = t_ppk_loc_pmin(1);
    end
end
ppk = MaxValue - MinValue;
timeFeature = [ppk, ppk_loc_pmax, ppk_loc_pmin];

% feature_vec = [feature_vec, timeFeature];
%% waveform features
[~, len] = size(data);
meanValue = mean(data, 2);tmp_mat = repmat(meanValue, [1, len]);
prepro_data = data - tmp_mat; % Y-direction shift to x-axis

waveFeature = extractWaveFeature(prepro_data);
feature_vec = [feature_vec, waveFeature];
%% feature normalization
max_feature = max(feature_vec, [], 1); min_feature = min(feature_vec, [], 1);% normalize
feature_vec = (feature_vec - repmat(min_feature,[N,1]))./repmat((max_feature - min_feature),[N,1]);
%feature_vec = feature_vec(N:-1:1,:); % Align with concentration vector

%% using regression algorithm to predict
% 'Evidential', 'Linear', 'NeuralNet', 'SVR', 'RegTree', 'LS-SVR'
select_regression_alg = {'Evidential','Linear', 'NeuralNet', 'SVR', 'RegTree', 'LS-SVR'};

[TotalNum,N_feature] = size(feature_vec);

% N_div=10;
% L=floor(N/(N_feature)); div_vec = 0.5:((1-0.5)/(N_div-1)):1; % = [0.25, 0.5, 0.75, 1.0];

k = 5;
foldSizes = computeFoldSizes(TotalNum, k);
disp("------------------------------------------------");
fprintf("[k = %d] k fold validation\n", k);
accuracy_vec = zeros(k, 1);
% this operation only for using the k-cross-validation code
[X_sorted, y_sorted] = randSortAndGroup([feature_vec, concentration_list'], ones(TotalNum, 1), 1);

fit_mean_mat = zeros(length(select_regression_alg), k); fit_std_mat = zeros(length(select_regression_alg), k);
fit_rmse_mat = zeros(length(select_regression_alg), k);
for roundNumber = 1:k
    y_sorted = [y_sorted,(1:N)'];  
    [X_train, tt, X_val, ~] = getFoldVectors(X_sorted, y_sorted, 1, TotalNum, foldSizes, roundNumber);
    torder = zeros(TotalNum, 1);torder(tt(:,2)) = 1;
    Y_train = X_train(:,N_feature+1); Y_val = X_val(:, N_feature+1);
    X_train = X_train(:,1:N_feature); X_val = X_val(:,1:N_feature);

    %% regression parts
    if sum(strcmpi(select_regression_alg, 'Evidential')) == 1 % using evidental regression
        [vec_evid] = RegressionTest_kcross(feature_vec, concentration_list, torder, 1, 'Evidential');
        fit_mean_evid = mean(vec_evid(1,:)); fit_std_evid = mean(vec_evid(2,:)); fit_rmse_evid = mean(vec_evid(3,:));        
        % for save
        Index_tmp = strcmpi(select_regression_alg, 'Evidential');
        fit_mean_mat(Index_tmp,roundNumber) = fit_mean_evid;
        fit_std_mat(Index_tmp,roundNumber) = fit_std_evid;
        fit_rmse_mat(Index_tmp,roundNumber) = fit_rmse_evid;
    end

    if sum(strcmpi(select_regression_alg, 'NeuralNet')) == 1 % using Nerual network regression
        [vec_net]  = RegressionTest_kcross(feature_vec, concentration_list, torder, 1, 'NeuralNet');
        fit_mean_net = mean(vec_net(1,:)); fit_std_net = mean(vec_net(2,:)); fit_rmse_net = mean(vec_net(3,:));        
        % for save
        Index_tmp = strcmpi(select_regression_alg, 'NeuralNet');
        fit_mean_mat(Index_tmp,roundNumber) = fit_mean_net;
        fit_std_mat(Index_tmp,roundNumber) = fit_std_net;
        fit_rmse_mat(Index_tmp,roundNumber) = fit_rmse_net;
    end

    if sum(strcmpi(select_regression_alg, 'Linear')) == 1 % using linear regression
        [vec_linear]  = RegressionTest_kcross(feature_vec, concentration_list, torder, 1, 'Linear');
        fit_mean_linear = mean(vec_linear(1,:)); fit_std_linear = mean(vec_linear(2,:)); fit_rmse_linear = mean(vec_linear(3,:));
        % for save
        Index_tmp = strcmpi(select_regression_alg, 'Linear');
        fit_mean_mat(Index_tmp,roundNumber) = fit_mean_linear;
        fit_std_mat(Index_tmp,roundNumber) = fit_std_linear;
        fit_rmse_mat(Index_tmp,roundNumber) = fit_rmse_linear;
    end

    if sum(strcmpi(select_regression_alg, 'SVR')) == 1 % using SVR
        [vec_svr]  = RegressionTest_kcross(feature_vec, concentration_list, torder, 1, 'SVR');
        fit_mean_svr = mean(vec_svr(1,:)); fit_std_svr = mean(vec_svr(2,:)); fit_rmse_svr = mean(vec_svr(3,:));
        % for save
        Index_tmp = strcmpi(select_regression_alg, 'SVR');
        fit_mean_mat(Index_tmp,roundNumber) = fit_mean_svr;
        fit_std_mat(Index_tmp,roundNumber) = fit_std_svr;
        fit_rmse_mat(Index_tmp,roundNumber) = fit_rmse_svr;
    end

    if sum(strcmpi(select_regression_alg, 'RegTree')) == 1 % using Regression Tree
        [vec_tree]  = RegressionTest_kcross(feature_vec, concentration_list, torder, 1, 'RegTree');
        fit_mean_tree = mean(vec_tree(1,:)); fit_std_tree = mean(vec_tree(2,:)); fit_rmse_tree = mean(vec_tree(3,:));       
        % for save
        Index_tmp = strcmpi(select_regression_alg, 'RegTree');
        fit_mean_mat(Index_tmp,roundNumber) = fit_mean_tree;
        fit_std_mat(Index_tmp,roundNumber) = fit_std_tree;
        fit_rmse_mat(Index_tmp,roundNumber) = fit_rmse_tree;      
    end

    if sum(strcmpi(select_regression_alg, 'LS-SVR')) == 1 % using LS-SVR
        [vec_lssvr]  = RegressionTest_kcross(feature_vec, concentration_list, torder, 1, 'LS-SVR');
        fit_mean_lssvr = mean(vec_lssvr(1,:)); fit_std_lssvr = mean(vec_lssvr(2,:)); fit_rmse_lssvr = mean(vec_lssvr(3,:));
        % for save
        Index_tmp = strcmpi(select_regression_alg, 'LS-SVR');
        fit_mean_mat(Index_tmp,roundNumber) = fit_mean_lssvr;
        fit_std_mat(Index_tmp,roundNumber) = fit_std_lssvr;
        fit_rmse_mat(Index_tmp,roundNumber) = fit_rmse_lssvr;      
    end

    if sum(strcmpi(select_regression_alg, 'CNN-net')) == 1 % using CNN-net
        [vec_cnn]  = RegressionTest_kcross(feature_vec, concentration_list, torder, 1, 'CNN-net'); % only repeat 1 time
        fit_mean_cnn = mean(vec_cnn(1,:)); fit_std_cnn = mean(vec_cnn(2,:)); fit_rmse_cnn = mean(vec_cnn(3,:));
        % for save
        Index_tmp = strcmpi(select_regression_alg, 'CNN-net');
        fit_mean_mat(Index_tmp,roundNumber) = fit_mean_cnn;
        fit_std_mat(Index_tmp,roundNumber) = fit_std_cnn;
        fit_rmse_mat(Index_tmp,roundNumber) = fit_rmse_cnn;
    end
end
%% show the result
if sum(strcmpi(select_regression_alg, 'Evidential')) == 1 % using evidental regression
    fprintf("[evid-predict]  :mean_value %f; std_value %f; RMSE %f\n", mean(fit_mean_evid), mean(fit_std_evid), mean(fit_rmse_evid));
end
if sum(strcmpi(select_regression_alg, 'NeuralNet')) == 1 % using Nerual network regression
    fprintf("[net-predict]   :mean_value %f; std_value %f; RMSE %f\n", mean(fit_mean_net), mean(fit_std_net), mean(fit_rmse_net));
end
if sum(strcmpi(select_regression_alg, 'Linear')) == 1 % using linear regression
    fprintf("[linear-predict]:mean_value %f; std_value %f; RMSE %f\n", mean(fit_mean_linear), mean(fit_std_linear), mean(fit_rmse_linear));        
end
if sum(strcmpi(select_regression_alg, 'SVR')) == 1 % using SVR
    fprintf("[SVR-predict]   :mean_value %f; std_value %f; RMSE %f\n", mean(fit_mean_svr), mean(fit_std_svr), mean(fit_rmse_svr));
end
if sum(strcmpi(select_regression_alg, 'RegTree')) == 1 % using Regression Tree
    fprintf("[tree-predict]  :mean_value %f; std_value %f; RMSE %f\n", mean(fit_mean_tree), mean(fit_std_tree), mean(fit_rmse_tree));
end
if sum(strcmpi(select_regression_alg, 'LS-SVR')) == 1 % using LS-SVR
    fprintf("[LSSVR-predict]  :mean_value %f; std_value %f; RMSE %f\n", mean(fit_mean_lssvr), mean(fit_std_lssvr), mean(fit_rmse_lssvr));
end
if sum(strcmpi(select_regression_alg, 'CNN-net')) == 1 % using CNN-net
    fprintf("[CNN-predict]  :mean_value %f; std_value %f; RMSE %f\n", mean(fit_mean_cnn), mean(fit_std_cnn), mean(fit_rmse_cnn));
end

