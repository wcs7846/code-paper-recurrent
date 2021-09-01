function [ fitTrain, fitTest, isConverged ] = svr( train_data, test_data, range )
%% SVR: support vector machine regression model
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  train_data  --> train data
         test_data   --> test data
         range       --> the select range (dim: 1D)
 Output: fitTrain    --> the fitting result of train data
         fitTest     --> the fitting result of test data
         isConverged --> a flag 

 Tips: train_data.data (dim:2D) the train vectors are row vector
       train_data.value(dim:1D) the responding value of train vectors
       the form of test_data is same as train_data
%}
[~, col] = size(train_data.data);
norm_cf = train_data.data./repmat(max(train_data.data,[],2), [1, col]);
Mdl = fitrsvm(norm_cf(:,range), train_data.value);
if Mdl.ConvergenceInfo.Converged
%     disp(sprintf('Suceess Converged'));
    isConverged = 1;
    %% test data
    norm_test_cf = test_data.data./repmat(max(test_data.data,[],2), [1, col]);
    fitTest = predict(Mdl, norm_test_cf(:,range));
    
    %% train data  
    norm_train_cf = train_data.data./repmat(max(train_data.data,[],2), [1, col]);
    fitTrain = predict(Mdl, norm_train_cf(:,range));
else
    isConverged = 0;
    fitTrain = [];
    fitTest = [];
%     disp(sprintf('Failed Converged'));
end

end

