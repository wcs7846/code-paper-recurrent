function [ output_args ] = RegressionTest_kcross( feature_vec, concentration_list, torder, K, AlgName, AlgPara )
%REGRESSIONTEST 此处显示有关此函数的摘要
%   此处显示详细说明
%   x：特征矩阵
%   y：回归目标
%   torder：分割训练集与测试集/ 1->训练(train); 0->测试(test);
%   K：重复次数
%   AlgName: 回归算法名称
%   AlgName: 回归算法参数
%   output_args: [误差均值 误差方差 RMSE] = 2*K
%---------------------------------------------------
% AlgName = 'Linear', 'NeuralNet', 'Evidential', 'SVR', 'RegTree', 'LS-SVR';
%% parameter verification
if nargin == 4 % for input 
    error('[Error]: Please enter the regression algorithm name!');
elseif nargin == 5
    % 'Linear' and 'SVR' and 'RegTree' is not need
    % 'NeuralNet'
    if strcmpi(AlgName, 'NeuralNet') || strcmpi(AlgName, 'NN')
        AlgPara.netpoint = 50;
    % 'Evidential'
    elseif strcmpi(AlgName, 'Evidential')
        AlgPara.alpha0 = 0.95;
        AlgPara.g0 = 0.5;
        AlgPara.w = ones(length(torder),1);
    elseif strcmpi(AlgName, 'LS-SVR')
        AlgPara.C = 100;
        AlgPara.g = 0.01;
    end
end

for k = 1:K
    N=length(torder);
    %% select the corresponding model and training
    if strcmpi(AlgName, 'Linear')
        mdl = fitlm(feature_vec(torder==1,:),concentration_list(torder==1)','linear','RobustOpts','on');
    elseif strcmpi(AlgName, 'NeuralNet') || strcmpi(AlgName, 'NN')
        net = fitnet(AlgPara.netpoint);
        net = train(net, feature_vec(torder==1,:)',concentration_list(torder==1));
    elseif strcmpi(AlgName, 'Evidential')
        alpha0=AlgPara.alpha0;  % input parameters for tbmfit
        g0=AlgPara.g0;
        ymin=min(concentration_list);  % min and max of output domain
        ymax=max(concentration_list);
        w=AlgPara.w;  % generation of weights (here, they are all taken equal to 1)
        
        [g,err] = tbmreg_fit(feature_vec(torder==1,:),concentration_list(torder==1)',w,N,alpha0,g0,ymin,ymax); % optimization of parameter g
    elseif strcmpi(AlgName, 'SVR')
        mdl = fitrsvm(feature_vec(torder==1,:),concentration_list(torder==1)','Standardize',true);
    elseif strcmpi(AlgName, 'RegTree') % Regression Tree
        [~,Nattr] = size(feature_vec);
        mdl = fitrtree(feature_vec(torder==1,:),concentration_list(torder==1)','MinParent',30);
    elseif strcmpi(AlgName, 'LS-SVR')
        % each sample is col vector
        [ alfa, bias ] = lssvr_train( feature_vec(torder==1,:)', concentration_list(torder==1)', AlgPara.C, AlgPara.g );
    elseif strcmpi(AlgName, 'CNN-net')
        layers = [ ...
            imageInputLayer([6 1 1])
            convolution2dLayer(12,25)
            reluLayer
            fullyConnectedLayer(1)
            regressionLayer];
        options = trainingOptions('sgdm','InitialLearnRate',0.001, ...
            'MaxEpochs',15);
        net = trainNetwork(feature_vec(torder==1,:)',concentration_list(torder==1)',layers,options);

    end
    %% predict
    if strcmpi(AlgName, 'Linear') || strcmpi(AlgName, 'SVR')
        predict_res = predict(mdl,feature_vec(torder==0,:));
    elseif strcmpi(AlgName, 'NeuralNet') || strcmpi(AlgName, 'NN')
        predict_res = net(feature_vec(torder==0,:)');
    elseif strcmpi(AlgName, 'Evidential')
        [mt,mtn,predict_res,z] = tbmreg_val(feature_vec(torder==1,:),concentration_list(torder==1)',w,feature_vec(torder==0,:), ...
            N,alpha0,g,ymin,ymax);  % computation of predictions for test examples
    elseif strcmpi(AlgName, 'RegTree')
        predict_res = predict(mdl,feature_vec(torder==0,:));
    elseif strcmpi(AlgName, 'LS-SVR')
        % each sample is col vector
        predict_res = lssvr_predict( feature_vec(torder==0,:)', feature_vec(torder==1,:)', alfa, bias, 0.01);
    elseif strcmpi(AlgName, 'CNN-net')
        predict_res = predict(net,feature_vec(torder==0,:)');
    end

    mean_value_res = mean(abs(predict_res - concentration_list(torder==0)'));
    var_value_res  = std(abs(predict_res - concentration_list(torder==0)'));
    % RMSE
    rmse = sqrt(sum((predict_res - concentration_list(torder==0)').^2)/N);
    %% output
    output_args(:,k) = [mean_value_res, var_value_res, rmse];
    
end

end

