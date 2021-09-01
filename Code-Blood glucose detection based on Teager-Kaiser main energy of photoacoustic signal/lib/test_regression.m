function [error_train, error_train_r, error_test, error_test_r] = test_regression( dataset, ...
          concentration_list, proprety_type, regression_type, range, debug_RegressionAnalysis)
%TEST_REGRESSION 此处显示有关此函数的摘要
%   此处显示详细说明
repeat_times = 1;
error_train_list = zeros(repeat_times, 1);
error_test_list  = zeros(repeat_times, 1);
for n = 1:repeat_times
    max_c = max(concentration_list); min_c = min(concentration_list);
    [trainset, testset] = divideDataset(dataset, concentration_list, 13);
    
    % use the SVR to fit the data    
    if strcmpi(proprety_type, 'tk') % Teager energy
        train_data.data = trainset.tk;
        test_data.data = testset.tk;
    elseif strcmpi(proprety_type, 'rmsf') % root means of frequency
        train_data.data = trainset.rmsf;
        test_data.data = testset.rmsf;
    elseif strcmpi(proprety_type, 'cf') % center frequency
        train_data.data = trainset.cf;
        test_data.data = testset.cf;
    elseif strcmpi(proprety_type, 'ibw') %
        train_data.data = trainset.ibw;
        test_data.data = testset.ibw;
    elseif strcmpi(proprety_type, 'fag')
        train_data.data = trainset.fag;
        test_data.data = testset.fag;
    elseif strcmpi(proprety_type, 'original') % original signal
        train_data.data = trainset.original;
        test_data.data = testset.original;
    else
        error('This proprety is not exisit!!');
    end
    
    train_data.value = trainset.concentraion;
    test_data.value = testset.concentraion;
    
    pp = range;
%     pp =  [1457, 1538];%pp = 1000:2000;
    if strcmpi(regression_type, 'linear')
        [fitTrain, fitTest, isConverged] = lrm(train_data, test_data, pp);
    elseif strcmpi(regression_type, 'svr')
        [fitTrain, fitTest, isConverged] = svr(train_data, test_data, pp);
    end
    %   [fitTrain, fitTest, isConverged] = svr(train_data, test_data, pp);
    
    if ~isConverged
        continue;
    end
    
    if debug_RegressionAnalysis && (n==1)
        if isConverged
            disp(sprintf('[%d]Suceess Converged', n));
        else
            disp(sprintf('[%d]Failed Converged', n));
        end
        
        figure;
        plot(1:length(testset.concentraion), testset.concentraion, 'ro-');hold on;
        plot(1:length(testset.concentraion), fitTest, 'bo-');hold on;
        legend('Real', 'Predict');title('test set comp');
        
        figure;
        plot(1:length(trainset.concentraion), trainset.concentraion, 'ro-');hold on;
        plot(1:length(trainset.concentraion), fitTrain, 'bo-');hold on;
        legend('Real', 'Predict');title('train set comp');
    end
    % calculate the error
    error_train = sum(abs(trainset.concentraion - fitTrain'))/length(trainset.concentraion);
    error_test  = sum(abs(testset.concentraion - fitTest'))/length(testset.concentraion);
    
    error_train_list(n) = error_train;
    error_test_list(n) = error_test;
    
    if debug_RegressionAnalysis
        error_train_r = (error_train *(max_c - min_c));
        error_test_r = (error_test *(max_c - min_c));
        disp(sprintf('[%d] the trainset error: %f[%f mg/L]; the testset error: %f[%f mg/L];', ...
            n, ...
            error_train, error_train_r, ...
            error_test,  error_test_r));
    end
    
%     order = [1:2:length(concentration_list), 2:2:length(concentration_list)];fit_vec = zeros(1, length(concentration_list));
%     size_train_set = 13; size_test_set  = length(concentration_list) - size_train_set;
%     fit_vec(order(1:size_train_set)) = fitTrain;
%     fit_vec(order(size_test_set+1:length(concentration_list))) = fitTest;
%     fit_vec = fit_vec*(max_c - min_c) + min(concentration_list);
%     save('fitting.mat', 'fit_vec', 'concentration_list' );
end
average_err_train = mean(error_train_list);
average_err_test  = mean(error_test_list);

error_train_r = (average_err_train *(max_c - min_c)) ;
error_test_r = (average_err_test *(max_c - min_c));

% disp(strcat(sprintf('SNB = %f', SNB_list(SNB_select)) ,' dB'));
% disp(sprintf('the average trainset error: %f[%f mg/L]; the average testset error: %f[%f mg/L];', ...
%     error_train, error_train_r, ...
%     error_test,  error_test_r));

end

