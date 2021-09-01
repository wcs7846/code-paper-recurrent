function [ output_args ] = KernelCalibration_predict( x_train, x_test, type, para, alpha )
%KERNELCALIBRATION_PREDICT 此处显示有关此函数的摘要
%   此处显示详细说明
% [Input]
% x_train: data matrix - [Nsample, Nfeature](train)
% y_train: reference result - column vector
% x_test:  data matrix - [Nsample, Nfeature](test)
% type: the type of kernel == {'poly', 'gauss'}
% para: type == 'poly', para = d(dimension) ; type == 'gauss', type = sita;
% alpha:  coefficient column vector
% [Output]
% output_args: predict result - row vector

% Ref: [1] Pai P P, De A, Banerjee S. Accuracy Enhancement for Noninvasive Glucose 
% Estimation Using Dual-Wavelength Photoacoustic Measurements and Kernel-Based Calibration[J].
% Ieee Transactions on Instrumentation and Measurement, 2018, 67(1): 126-36
[Nsample_test, ~] = size(x_test); [Nsample_train,~] = size(x_train);
if strcmpi(type, 'poly') % using Polynomial kernel
    K_mat = (1 + x_test*x_train').^para; % poly
elseif strcmpi(type, 'gauss') % using Gaussian kernel
    tmp_K = abs(repmat(x_test, [Nsample_train, 1]) - repelem(x_train, Nsample_test, 1)); % ||pi - pj||_2
    tmp_K = sum(tmp_K.^2, 2);
    K_mat = exp(-reshape(tmp_K, [Nsample_train, Nsample_test])'./(2*para));
end
output_args = (K_mat*alpha)';

end

