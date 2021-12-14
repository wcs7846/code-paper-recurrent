function [ output_args ] = evidental_regression( feature_vec, concentration_list, torder, K )
%EVIDENTAL_REGRESSION 此处显示有关此函数的摘要
%   此处显示详细说明
%   x：特征矩阵
%   y：回归目标
%   torder：分割训练集与测试集
%   K：重复次数
%   output_args: [误差均值 误差方差 RMSE] = 2*K
for k = 1:K
    alpha0=0.95;  % input parameters for tbmfit
    g0=0.5;
    ymin=min(concentration_list);  % min and max of output domain
    ymax=max(concentration_list);
    w=ones(length(torder),1);  % generation of weights (here, they are all taken equal to 1)
    N=length(torder);
    
    [g,err] = tbmreg_fit(feature_vec(torder==1,:),concentration_list(torder==1)',w,N,alpha0,g0,ymin,ymax); % optimization of parameter g
    if length(torder) == sum(torder)
        [mt,mtn,predict_evid,z] = tbmreg_val(feature_vec(torder==1,:),concentration_list(torder==1)',w,feature_vec(torder==1,:),N,alpha0,g,ymin,ymax);  % computation of predictions for test examples
        mean_value_evid = mean(abs(predict_evid - concentration_list(torder==1)'));
        var_value_evid  = std(abs(predict_evid - concentration_list(torder==1)'));
        % RMSE
        rmse = sqrt(sum((predict_evid - concentration_list(torder==1)').^2)/N);
    else
        [mt,mtn,predict_evid,z] = tbmreg_val(feature_vec(torder==1,:),concentration_list(torder==1)',w,feature_vec,N,alpha0,g,ymin,ymax);  % computation of predictions for test examples
        mean_value_evid = mean(abs(predict_evid - concentration_list'));
        var_value_evid  = std(abs(predict_evid - concentration_list'));
        % RMSE
        rmse = sqrt(sum((predict_evid - concentration_list').^2)/N);
    end

    output_args(:,k) = [mean_value_evid, var_value_evid, rmse];
end

end

