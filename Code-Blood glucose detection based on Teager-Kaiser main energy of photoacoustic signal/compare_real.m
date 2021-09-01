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
load('real-data');
concentration_list = tmp_con';

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
sensor_data = data;
for n = length(concentration_list):-1:1
    data = sensor_data(n,:);
    sample_rate = 50e6; % 50MHz
    
    xrange = (1:numel(data))*1/sample_rate;
    if debug_PropretyAnalysis
        figure(show_rd);plot(xrange, data, 'color', lineColor(n, :));hold on;
        title(strcat('the real blood glucose signal: ',num2str(concentration_list(n)),'mg/L'));
        xlabel('time/s');ylabel('pressure');
        
        label = [label, string(strcat(num2str(concentration_list(n)), ' mg/L'))];
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
%     [tf_Stran,ts,fs] = tfrwv(data',1:Ndata, Ndata);
%     fs = fs';
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
%     title(sprintf('normal S transform %d mg/L', concentration_list(n)));
    
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
        
        label = [label, string(strcat(num2str(concentration_list(n)), ' mg/L'))]; % string function must be used on MATLAB2016b
        
        pause(0.001);
    end
    
    % show the processing status
    disp(sprintf(strcat(num2str(concentration_list(n)), ' mg/L', ' is completed!')));
    
    special_loc = 3936; % unit: um
    loc = find(abs(xrange-special_loc) == min(abs(xrange-special_loc)));
    value = tk(loc);
end