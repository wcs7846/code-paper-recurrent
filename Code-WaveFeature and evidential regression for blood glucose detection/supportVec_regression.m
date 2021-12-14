function [ output_args ] = supportVec_regression( feature_vec, concentration_list, torder, K  )
%LINEAR_REGRESSION 此处显示有关此函数的摘要
%   此处显示详细说明
%   x：特征矩阵
%   y：回归目标
%   torder：分割训练集与测试集
%   K：重复次数
%   output_args: [误差均值 误差方差 RMSE] = 2*K
for k = 1:K
        N=length(torder);
    mdl = fitrsvm(feature_vec(torder==1,:),concentration_list(torder==1)','Standardize',true);
    if length(torder) == sum(torder)
        predict_svr = predict(mdl,feature_vec(torder==1,:));
        mean_value_svr = mean(abs(predict_svr - concentration_list(torder==1)'));
        var_value_svr  = std(abs(predict_svr - concentration_list(torder==1)'));
                % RMSE
        rmse = sqrt(sum((predict_svr - concentration_list(torder==1)').^2)/N);
    else
        predict_svr = predict(mdl,feature_vec);
        mean_value_svr = mean(abs(predict_svr - concentration_list'));
        var_value_svr  = std(abs(predict_svr - concentration_list'));
                % RMSE
        rmse = sqrt(sum((predict_svr - concentration_list').^2)/N);
    end
    output_args(:,k) = [mean_value_svr, var_value_svr, rmse];
end


end

