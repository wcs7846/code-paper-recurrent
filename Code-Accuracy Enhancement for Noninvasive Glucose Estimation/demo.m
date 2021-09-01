%% Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
% This script is designed to model the signal model of Blood glucose signal
% Author: MarkLHF
% Date: 2021-9-1
% Reference:

close all;
clear; format long;

[x,t] = simplefit_dataset;
% train
d = 3; sita = 100;
max_t = max(t);min_t = min(t);
t = (t - mean(t))./max_t;
[ alpha ] = KernelCalibration_train( x', t, 'poly', d );
% test
y = KernelCalibration_predict(x', x', 'poly', d, alpha);
y = (y - mean(y))./max(y);

figure;
plot(x,t,'ro');hold on;
plot(x,y,'b-*');hold on;
legend('test-data', 'fitting');