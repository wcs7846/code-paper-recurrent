%% setting
close all;
clear; format long;

addpath(genpath('D:\Matlab_program\Matlab_toolbox\tftb-0.2\mfiles'));
%% Copyright(C),UESTC, School of Information and Communication Engineering, IDIP

% switch
debug_PropretyAnalysis = 1;
debug_RegressionAnalysis = 1;
debug_ConcentrationAnalysis = 1;

addpath('./lib');
addpath(genpath('./evaluation'));
% load the sensor data
% D:/龙鸿峰/研究室 — — 血糖项目/simulate/simulate_result
path = 'D:/龙鸿峰/研究室 — — 血糖项目/simulate/simulate_result/acoustic/';
concentration_list = 500:100:3000; concentration_list=concentration_list/10; % mg/dL
SNB_list = [1:10:90]-1;
% SNB_list = [ 1,  3,  5,  8, 10,...
%             15, 20, 30, 50,]; % unit: dB
SNB_select = 1; 
%SNB_value = SNB_list(SNB_select);
SNB_value = 15;

dataset.cf  = [];
dataset.ibw = [];
dataset.tk = [];
dataset.rmsf = [];
dataset.fag = [];
dataset.original = [];

lineColor = linspecer(length(concentration_list)); label=[];

% figure allocation
if debug_PropretyAnalysis
    show_rd = figure; % real data figure;
%     show_feature1 = figure; % feature figure;
%     show_feature2 = figure;
    show_feature3 = figure;
%     show_feature4 = figure;
    n = 26;
end

for n = length(concentration_list):-1:1
    load(strcat(path, 'acoustic_skin_', num2str(concentration_list(n)*10), '.mat'));
    
    [row, col] = size(sensor_data);
    data = sensor_data(round(row/2),:);
    dt = kgrid.dt; % the interval time == the reciprocal
%     average_v = 1200.0 * 1e6; % um/s
%     d = dt*average_v; % 0.01 mm = 10 um

    %% down-sample -- using resample function
    set_samplerate = 100*1e6;
    t_total = dt*length(data); num_point_ds = round(t_total*set_samplerate);
    ds_data = resample(data, num_point_ds, length(data));
    % update
    data = ds_data; dt = 1/set_samplerate;
    
    %% Band-pass filter
    tt_data = fftshift(fft(data, 2*length(data)));
    low_pass = 3.5e6; high_pass = 10.5e6;
    tl = round(low_pass/(set_samplerate/2)*length(data));
    th = round(high_pass/(set_samplerate/2)*length(data));
    tt_data_left  = (length(data) - th + 1):(length(data) - tl + 1);
    tt_data_right = (length(data) + tl - 1):(length(data) + th - 1);
    tt_data_2 = zeros(1,2*length(data));
    tt_data_2(tt_data_left) = tt_data(tt_data_left);
    tt_data_2(tt_data_right) = tt_data(tt_data_right);
    data2 = real(ifft(ifftshift(tt_data_2)));
    data2 = data2(1:length(data));
    
    data = data2;
    
    thick_vec = [100, 1400, 2000, 400, 200, 1400];
    speed_vec = [1645, 1595, 1450, 1547, 1500, 1547];
    
%     xrange = kgrid.t_array; x_axis_label = 'time/s';
    [xrange ,data] = timedepthConvert(data, speed_vec, thick_vec, dt);
    x_axis_label = 'depth/um';

    %% synthetic signal
    power_signal = mean(data.^2); % average power
    power_noise = power_signal*10^(-SNB_value/10);
    power_noise_t = 10*log10(power_noise); % transform unit to dBw
    noise = wgn(1, length(data), power_noise_t);
    
    data = data + noise;
    
    if debug_PropretyAnalysis
        figure(show_rd);plot(xrange, data, 'color', lineColor(n, :));hold on;
        title(strcat('the simulate blood glucose signal: ',num2str(concentration_list(n)),'mg/dL'));
        xlabel(x_axis_label);ylabel('pressure');xlim([min(xrange),max(xrange)]);
    end
    %% the signal model
    
    addpath('./lib');
    
    check_freq = 44e6; % 9.0 MHz
    check_time = 1461;
    Ndata = length(data);
    % padding = 10; len = length(data);
    % data_pad = padarray(data, [0, padding], 'symmetric');
    % [Part 1] calculate the time-frequence spectrum
%     dt = kgrid.dt;
    sample_rate = 1/dt;
    %% S transform
    [tf_Stran, ts, fs] = st(data, 1);
    %% STFT 
%     [tf_Stran,ts,fs] = tfrstft(data',1:Ndata, Ndata);
%     loc = find(fs < 0); fs(loc) = []; tf_Stran(loc,:) = [];
%     fs = fs';
    %% Gabor transform
%     [tf_Stran,ts,fs] = tfrgabor([data,0]', (Ndata+1)/2, (Ndata+1)/2);
%     tf_Stran = tf_Stran(1:Ndata,:); fs = fs(1:Ndata);
%     fs = fs';
    %% WVD distribution
    [tf_Stran,ts,fs] = tfrwv(data',1:Ndata, Ndata);
    fs = fs';
    % frequency slice
    freq_vect = fs*sample_rate;
%         for m = 1:length(check_freq)
%             loc = find(abs(freq_vect-check_freq)==min(abs(freq_vect-check_freq)));
%             freq_slice = tf_Stran(loc, :);
%             figure(show_feature1);
%             plot(xrange, abs(freq_slice), 'color', lineColor(n, :)); hold on;
%             xlabel(x_axis_label);ylabel('time-frequency amplitude');xlim([min(xrange),max(xrange)]);
%             title(strcat('check frequency: ', num2str(check_freq/1e6), 'MHz'));
%             ylim([0, 0.85])
%         end

%     upper_freq = 1/dt/2;low_freq = 0;
    upper_freq = 13e6; % 13MHz
    low_freq = 2e6; % 2MHz
    freq_band = find((freq_vect <= upper_freq) & (freq_vect >= low_freq));
    tf_Stran_band = tf_Stran(freq_band, :);
    
    % debug for time-frequence transfrom
%     figure;imagesc(ts(freq_band)/sample_rate, fs(freq_band)*sample_rate, abs(tf_Stran_band));
%     set(gca,'ydir','normal');colorbar;xlabel('Time');ylabel('Frequence');
%     title(sprintf('normal S transform %d mg/dL', concentration_list(n)));
    
    % [part 2] extract the property analysis from time-frequence spectrum
    [ cf, rmsf, ibw ] = centerFreq(tf_Stran_band, fs(freq_band)*sample_rate);
%     cf = cf/(max(cf)+eps); rmsf = rmsf/(max(rmsf)+eps); ibw = ibw/(max(ibw)+eps);
    dataset.original = [dataset.original; data];
    dataset.cf = [dataset.cf; cf];
    dataset.rmsf = [dataset.rmsf; rmsf];
    dataset.ibw = [dataset.ibw; ibw];
    
    tk = teagerEnergy(tf_Stran_band);
    dataset.tk = [dataset.tk; tk];
    
    fag = freqAttenGrad(tf_Stran_band, freq_vect(freq_band)); 
    dataset.fag = [dataset.fag; abs(fag)];
    %     [cf, max_cf, min_cf] = normalization(cf);
    %     [rmsf, max_rmsf, min_rmsf] = normalization(rmsf);
    %     [ibw,  max_ibw,  min_ibw ] = normalization(ibw);
    
    %% show
    show_x = find(xrange > 000 & xrange < 5000 == 1) ;
    if debug_PropretyAnalysis
%         figure(show_feature1); plot(xrange(show_x), rmsf(show_x), 'color', lineColor(n, :)); hold on;
%         figure(show_feature2); plot(xrange(show_x), ibw(show_x),  'color', lineColor(n, :)); hold on;
        figure(show_feature3); plot(xrange(show_x), tk(show_x),   'color', lineColor(n, :)); hold on;
%         figure(show_feature4); plot(xrange(show_x), fag(show_x),  'color', lineColor(n, :)); hold on;
        
        %     [xt, yt] = meshgrid(ts/sample_rate,  fs*sample_rate);
        %     mesh(xt, yt, abs(tf_Stran)); hold on;
        %     set(gca,'ydir','normal');colorbar;xlabel('Time');ylabel('Frequence');
        
        label = [label, string(strcat(num2str(concentration_list(n)), ' mg/dL'))]; % string function must be used on MATLAB2016b
        xlabel(x_axis_label);ylabel('TKME');xlim([min(xrange),max(xrange)]);
        pause(0.001);
    end
    
    % show the processing status
    disp(sprintf(strcat(num2str(concentration_list(n)), ' mg/dL', ' is completed!')));
    
    special_loc = 3936; % unit: um
    loc = find(abs(xrange-special_loc) == min(abs(xrange-special_loc)));
    value = tk(loc);
end
xlim([min(xrange(show_x)),max(xrange(show_x))]);
colorbar('Ticks',0:1/(length(concentration_list)-1):1,'TickLabels', fliplr(label));colormap(lineColor);