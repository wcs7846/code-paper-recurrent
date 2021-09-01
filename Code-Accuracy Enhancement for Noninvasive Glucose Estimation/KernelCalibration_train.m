function [ alpha ] = KernelCalibration_train( x, y, type, para )
%KERNELCALIBRATION 此处显示有关此函数的摘要
%   此处显示详细说明
% [Input]
% x: data matrix - [Nsample, Nfeature]
% y: reference result - column vector
% type: the type of kernel == {'poly', 'gauss'}
% para: type == 'poly', para = d(dimension) ; type == 'gauss', type = sita;
% [Output]
% alpha: model coefficient, column vector

% Ref: [1] Pai P P, De A, Banerjee S. Accuracy Enhancement for Noninvasive Glucose 
% Estimation Using Dual-Wavelength Photoacoustic Measurements and Kernel-Based Calibration[J].
% Ieee Transactions on Instrumentation and Measurement, 2018, 67(1): 126-36
[Nsample, Nfeature] = size(x);
if strcmpi(type, 'poly') % using Polynomial kernel
    K = (1+x*x').^para;
elseif strcmpi(type, 'gauss') % using Gaussian kernel
    tmp_K = abs(repmat(x, [Nsample, 1]) - repelem(x, Nsample, 1)); % ||pi - pj||_2
    tmp_K = sum(tmp_K.^2, 2);
    K = exp(-reshape(tmp_K, [Nsample, Nsample])./(2*para));
end


alpha_0 = ones(Nsample, 1);

options = optimoptions(@fminunc,'Algorithm','quasi-newton');
lambda = 0.1;
fun = @(alpha)loss(alpha, K, y, lambda);
[alpha_opti, fval] = fminunc(fun,alpha_0,options);

alpha = alpha_opti;
end

function [f] = loss(alpha, K, y, lambda)
% Calculate objective f
f = norm((K*alpha - y),2) + lambda*alpha'*K*alpha;
% g = 2*K'*K*alpha-2*K'*y+lambda*(K'+K)*alpha;
end