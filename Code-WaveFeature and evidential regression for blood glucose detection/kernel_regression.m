function [ output_args ] = kernel_regression( feature_vec, concentration_list, torder, K )
%KERNEL_REGRESSION 此处显示有关此函数的摘要
%   此处显示详细说明
for k = 1:K
    % train
    d = 3; sita = 100;
    [ alpha ] = KernelCalibration_train( feature_vec(torder==1,:), concentration_list(torder==1)','poly', d );
    
    if length(torder) == sum(torder)
        predict_res = KernelCalibration_predict(feature_vec(torder==1,:), feature_vec(torder==1,:), 'poly', d, alpha);
        
        mean_value_kernel = mean(abs(predict_res - concentration_list(torder==1)'));
        var_value_kernel  = std(abs(predict_res - concentration_list(torder==1)'));
    else
        predict_res = KernelCalibration_predict(feature_vec(torder==1,:), feature_vec, 'poly', d, alpha);
        
        mean_value_kernel = mean(abs(predict_res - concentration_list'));
        var_value_kernel  = std(abs(predict_res - concentration_list'));
    end
    output_args(:,k) = [mean_value_kernel, var_value_kernel];
end

end

