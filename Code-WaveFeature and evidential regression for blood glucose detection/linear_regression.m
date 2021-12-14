function [ output_args ] = linear_regression( feature_vec, concentration_list, torder, K  )
%LINEAR_REGRESSION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%   x����������
%   y���ع�Ŀ��
%   torder���ָ�ѵ��������Լ�
%   K���ظ�����
%   output_args: [����ֵ ���� RMSE] = 2*K
for k = 1:K
    N=length(torder);
    mdl = fitlm(feature_vec(torder==1,:),concentration_list(torder==1)','linear','RobustOpts','on');
    if length(torder) == sum(torder)
        predict_linear = predict(mdl,feature_vec(torder==1,:));
        mean_value_linear = mean(abs(predict_linear - concentration_list(torder==1)'));
        var_value_linear  = std(abs(predict_linear - concentration_list(torder==1)'));
                % RMSE
        rmse = sqrt(sum((predict_linear - concentration_list(torder==1)').^2)/N);
    else
        predict_linear = predict(mdl,feature_vec);
        mean_value_linear = mean(abs(predict_linear - concentration_list'));
        var_value_linear  = std(abs(predict_linear - concentration_list'));
                % RMSE
        rmse = sqrt(sum((predict_linear - concentration_list').^2)/N);
    end
    output_args(:,k) = [mean_value_linear, var_value_linear, rmse];
end


end

