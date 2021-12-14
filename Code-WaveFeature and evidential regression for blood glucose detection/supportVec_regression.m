function [ output_args ] = supportVec_regression( feature_vec, concentration_list, torder, K  )
%LINEAR_REGRESSION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%   x����������
%   y���ع�Ŀ��
%   torder���ָ�ѵ��������Լ�
%   K���ظ�����
%   output_args: [����ֵ ���� RMSE] = 2*K
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

