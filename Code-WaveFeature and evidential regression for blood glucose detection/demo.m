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

% addpath('./lib');
% addpath(genpath('./evaluation'));
% load the sensor data
load('real-data');
concentration_list = tmp_con'/10;

N = length(concentration_list);
lineColor = linspecer(numel(concentration_list)); label=[];
lineColor = lineColor(N:-1:1,:);
% figure allocation
if debug_PropretyAnalysis
    show_rd = figure; % real data figure;
end
feature_vec= [];
%% time-frequency features extraction
tk_mat = []; rmsf_mat = []; ibw_mat = []; fag_mat = [];
for n = 1:N % each sample
    tmpdata = data(n,:);
    sample_rate = 50e6; % 50MHz
    
    %% synthetic signal
    SNB_value = 5;
    power_signal = mean(tmpdata.^2); % average power
    power_noise = power_signal*10^(-SNB_value/10);
    power_noise_t = 10*log10(power_noise); % transform unit to dBw
    noise = wgn(1, length(tmpdata), power_noise_t);
    
%     data(n,:) = tmpdata + noise;
    
    xrange = (1:numel(data(n,:)))*1/sample_rate;
    if debug_PropretyAnalysis
        figure(show_rd);plot(xrange, data(n,:), 'color', lineColor(n, :));hold on;
        title(strcat('the real blood glucose signal: ',num2str(concentration_list(n)),'mg/dL'));
        xlabel('time/s');ylabel('pressure');
        
        label = [label, string(strcat(num2str(concentration_list(n)), ' mg/dL'))];
    end
    
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
        figure(show_feature3); plot(xrange(show_x), tk(show_x),   'color', lineColor(n, :)); hold on;
    end
end
% time-frequency feature
special_loc = 2.76e-6; % unit: um
loc = find(abs(xrange-special_loc) == min(abs(xrange-special_loc)));
% 使用时频特征作为特征向量
tf_vec = [tk_mat(:, loc), rmsf_mat(:, loc), ibw_mat(:, loc)]; %
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
waveFeature = extractWaveFeature(data);
feature_vec = [feature_vec, waveFeature];
%%
% select tk as the data
% data_select = tk_mat(:, loc); % data_select = x; concentration_list = y
% data_select = data_select(length(data_select):-1:1);
% % normalize - tk attributes
% max_attr = max(data_select); min_attr = min(data_select);
% % data_select = (data_select - min_attr)/(max_attr - min_attr);
% % normalize - concentration
% max_c = max(concentration_list); min_c = min(concentration_list);
% % c = (concentration_list - min_c)/(max_c - min_c);
% c = concentration_list;
%
% % Mdl = fitlm(data_select, c);
% % pred_vec = predict(Mdl, data_select);
%
% data_ex = [data_select'; ones(1, length(data_select))]';
% para_vec = inv(data_ex'*data_ex)*data_ex'*c';
% Num_plot = 100;tt = min_attr-0.2*min_attr:(max_attr-min_attr)/(Num_plot-1):max_attr+0.2*max_attr;
% line_vec = zeros(1, numel(tt));
% for n = 1:numel(tt)
%     tmp = [tt(n), 1]*para_vec;
%     line_vec(n) = (tmp);
% end
% plot(tt, line_vec, 'b-');hold on;xlim([min(tt), max(tt)]);
% legend('real data','fitting curve');
%
% pred_vec = zeros(1, length(data_select));
% for n = 1:length(data_select)
%     tmp = [data_select(n), 1]*para_vec;
%     pred_vec(n) = tmp;
% end
%
% % show
% figure;
% plot(1:length(concentration_list), concentration_list, 'r+');hold on;
% plot(1:length(concentration_list), pred_vec, 'bo'); hold on;
% legend('real', 'predict');
% save('fitting_data.mat', 'pred_vec', 'concentration_list');
%
% %
% mean_value = mean(pred_vec - concentration_list);
% var_value  = std(pred_vec - concentration_list);
% re = (pred_vec - concentration_list)./concentration_list;
%
% mean_value_re = mean(re);
% var_value_re = std(re);
[~,N_feature] = size(feature_vec);N_div=10;
L=floor(N/(N_feature+2)); div_vec = 1/L:((1-1/L)/(N_div-1)):1; % = [0.25, 0.5, 0.75, 1.0];
fit_mean_mat = zeros(4, N_div); fit_std_mat = zeros(4, N_div);
show_Reg_mean = figure;     show_Reg_var = figure;
for l = 1:N_div
    %% using fitnet to fitting the result
    disp(sprintf("---[Prop:%f]---", div_vec(l)*100));
    % Repeat the experiment K times
    K = 10;
    vec_net = zeros(2, K); % [mean_value; std_value] for one time
    vec_evid = zeros(2, K);
    vec_linear = zeros(2, K);
    vec_kernel = zeros(2, K);
    if div_vec(l) == 1
        torder = [ones(1, N)]; % select all as train
    else
        torder = [ones(1, round(N*div_vec(l))), zeros(1, N-round(N*div_vec(l)))]; % select 1/l as test
    end
    randIndex = randperm(size(torder, 2));
    torder = torder(randIndex);
    
    % 使用整个数据作为特征向量
    % net = fitnet(100);
    % net = train(net, data(torder==1,:)',concentration_list(torder==1));
    % predict_net = net(data(torder==0,:)');
    
    % feature_vec = data;
    max_feature = max(feature_vec, [], 1); min_feature = min(feature_vec, [], 1);% normalize
    feature_vec = (feature_vec - repmat(min_feature,[N,1]))./repmat((max_feature - min_feature),[N,1]);
    %feature_vec = feature_vec(N:-1:1,:); % Align with concentration vector
    
    %% using evidental regression
    [ vec_evid ] = evidental_regression( feature_vec, concentration_list, torder, K );
    %% using Nerual network regression
    [ vec_net ] = nn_regression( feature_vec, concentration_list, torder, K );    
    %% using linear regression
    [ vec_linear ] = linear_regression( feature_vec, concentration_list, torder, K  );
    %% using SVR
    [ vec_svr] = supportVec_regression( feature_vec, concentration_list, torder, K);
    %% using kernel calibration
%     [ vec_kernel ] = kernel_regression( feature_vec, concentration_list, torder, K );
    
    fit_mean_mat(:,l)= [mean(vec_net(1,:)); mean(vec_evid(1,:)); mean(vec_linear(1,:)); mean(vec_svr(1,:))];
    fit_std_mat(:,l) = [mean(vec_net(2,:)); mean(vec_evid(2,:)); mean(vec_linear(2,:)); mean(vec_svr(2,:))];
    fit_rmse_mat(:,l) = [mean(vec_net(3,:)); mean(vec_evid(3,:)); mean(vec_linear(3,:)); mean(vec_svr(3,:))];
    
    disp(sprintf("[linear-predict]:mean_value %f; std_value %f; RMSE %f", mean(vec_linear(1,:)), mean(vec_linear(2,:)), mean(vec_linear(3,:))));
    disp(sprintf("[SVR-predict]:mean_value %f; std_value %f; RMSE %f", mean(vec_svr(1,:)), mean(vec_svr(2,:)), mean(vec_svr(3,:))));
    disp(sprintf("[net-predict]   :mean_value %f; std_value %f; RMSE %f", mean(vec_net(1,:)), mean(vec_net(2,:)), mean(vec_net(3,:))));
    disp(sprintf("[evid-predict]  :mean_value %f; std_value %f; RMSE %f", mean(vec_evid(1,:)), mean(vec_evid(2,:)), mean(vec_evid(3,:))));
%     disp(sprintf("[kernel-predict]:mean_value %f; std_value %f; RMSE %f", mean(vec_kernel(1,:)), mean(vec_kernel(2,:)), mean(vec_net(3,:))));
end
%% draw figure
figure(show_Reg_mean);
setLineWidth = 2;
plot(div_vec*100,fit_mean_mat(1,:),'r-*', 'LineWidth', setLineWidth);hold on;
plot(div_vec*100,fit_mean_mat(3,:),'k-*', 'LineWidth', setLineWidth);hold on;
% plot(div_vec*100,fit_mean_mat(4,:),'g-*', 'LineWidth', setLineWidth);hold on;
plot(div_vec*100,fit_mean_mat(4,:),'g-*', 'LineWidth', setLineWidth);hold on;
plot(div_vec*100,fit_mean_mat(2,:),'b-*', 'LineWidth', setLineWidth);hold on;
legend('Neural Network','linear','SVR','Our method');
xlabel('train proportion(%)');ylabel('Mean value of error(mg/dL)');
xlim([div_vec(1)*100, 100]);

figure(show_Reg_var);
plot(div_vec*100,fit_std_mat(1,:),'r-*', 'LineWidth', setLineWidth);hold on;
plot(div_vec*100,fit_std_mat(3,:),'k-*', 'LineWidth', setLineWidth);hold on;
% plot(div_vec*100,fit_std_mat(4,:),'g-*', 'LineWidth', setLineWidth);hold on;
plot(div_vec*100,fit_std_mat(4,:),'g-*', 'LineWidth', setLineWidth);hold on;
plot(div_vec*100,fit_std_mat(2,:),'b-*', 'LineWidth', setLineWidth);hold on;
legend('Neural Network','linear','SVR','Our method');
xlabel('train proportion(%)');ylabel('Variance of error(mg/dL)');
xlim([div_vec(1)*100, 100]);