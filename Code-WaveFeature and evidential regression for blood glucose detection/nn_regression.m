function [ output_args ] = nn_regression( feature_vec, concentration_list, torder, K )
%NN_REGRESSION �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%   x����������
%   y���ع�Ŀ��
%   torder���ָ�ѵ��������Լ�
%   K���ظ�����
%   output_args: [����ֵ ���� RMSE] = 2*K
for k = 1:K
    net = fitnet(50);    N=length(torder);
    net = train(net, feature_vec(torder==1,:)',concentration_list(torder==1)); %
    if length(torder) == sum(torder)
        predict_net = net(feature_vec(torder==1,:)');
        mean_value_net = mean(abs(predict_net - concentration_list(torder==1)));
        var_value_net  = std(abs(predict_net - concentration_list(torder==1)));
                % RMSE
        rmse = sqrt(sum((predict_net - concentration_list(torder==1)').^2)/N);
    else
        predict_net = net(feature_vec');
        mean_value_net = mean(abs(predict_net - concentration_list));
        var_value_net  = std(abs(predict_net - concentration_list));
                % RMSE
        rmse = sqrt(sum((predict_net - concentration_list').^2)/N);
    end
    output_args(:,k) = [mean_value_net, var_value_net, rmse];
end

end

