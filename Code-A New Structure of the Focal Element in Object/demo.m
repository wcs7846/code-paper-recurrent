%% Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
% This script is designed to test the information fractals for evidential
% pattern classification
% Author: MarkLHF
% Date: 2020-1-7(first version)
%       2020-10-20
%       2020-11-17
%% Experiment 1: Constructive database test
% traditional focal element + traditional classification method
clc;clear all;close all;

addpath('./lib')
debug_showdistribution = 0;
debug_kfoldValidation = 1;
%% use parpool to quicken the calculation
% parpool;
% 'iris'; 'wine'; 'wdbc'; 'sonar'; 'ionosphere'
select = 'iris';
%% IRIS dataset
if strcmpi(select, 'iris')
    %% load the iris.csv
    disp("------------------------------------------------");
    disp("Using the IRIS dataset to test");
    fileID = fopen('iris_matlab.txt');
    data = textscan(fileID,'%f\t%f\t%f\t%f\t%f\t%s');
    fclose(fileID);
    
    if ~isempty(data)
        N = length(data{1});
    end
    typeClass = [{'setosa'}; {'versicolor'}; {'virginica'}];
    typeAttr  = [{'Sepal_len'}; {'Sepal_wid'}; {'Petal_len'}; {'Petal_wid'}];
    Num_class = length(typeClass);
    Num_Attr  = length(typeAttr);
    % construct the database
    tdata = [ data{2}, data{3}, data{4}, data{5} ]; ttype = zeros(N,1);
    TotalNum = zeros(Num_class, 1);
    for n = 1:N
        for m = 1:numel(typeClass)
            if strcmpi(data{6}{n}, typeClass{m})
                ttype(n) = m;
                break;
            end            
        end
        TotalNum(ttype(n)) = TotalNum(ttype(n)) + 1;
    end
    data = [ttype, tdata];
    %% transform the source data to BPA
    evid_set = loadIRIS_new(data, TotalNum, typeClass);
end
%% WINE dataset
if strcmpi(select, 'wine')
    %% load the wine.data
    disp("------------------------------------------------");
    disp("Using the WINE dataset to test");
    
    data=load('wine.data');
    if ~isempty(data)
        [N, tt] = size(data);
    end   
    Num_Attr  = tt-1; % The first is class;
    
    typeClass = [{'c1'}; {'c2'}; {'c3'}];
    
    Num_class = length(typeClass);
    [ evid_set ] = loadWine_new( data, typeClass );
    TotalNum = zeros(Num_class, 1);
    for n = 1:N
        tmp = data(n,1);
        TotalNum(tmp) = TotalNum(tmp) + 1;
    end
end
%% Wisconsin Diagnostic Breast Cancer (WDBC) dataset
if strcmpi(select, 'wdbc')
    %% load the wdbc.data
    disp("------------------------------------------------");
    disp("Using the WDBC dataset to test");

    data = importfile_wdbc('wdbc.data', 1, 569);
    % simple transform: M --> 1; B --> 2;
    if ~isempty(data)
        [N, tt] = size(data);
    end 
    tmp_data = zeros(N, tt - 1); ID_vec = zeros(N,1);
    for n = 1:N
        ID_vec(n) = data{n,1};
        tclass = data{n,2};
        if strcmpi(tclass, 'M')
            tmp_data(n,1) = 1;
        elseif strcmpi(tclass, 'B')
            tmp_data(n,1) = 2;
        end
        % other data
        for m = 3:tt
            tmp_data(n,m-1) = data{n,m};
        end        
    end
       
    % re-order
    s1 = tmp_data(tmp_data(:,1) == 1, :);
    s2 = tmp_data(tmp_data(:,1) == 2, :);
    data = [s1;s2];
%     data = tmp_data;
    if ~isempty(data) %^update
        [N, tt] = size(data);
    end
    
    Num_Attr  = tt-1; % The first is class;
    
    typeClass = [{'malignant'}; {'benign'}];
    Num_class = length(typeClass);
    [ evid_set ] = loadWDBC_new( data, typeClass );
    

    
    TotalNum = zeros(Num_class, 1);
    for n = 1:N
        tmp = data(n,1);
        TotalNum(tmp) = TotalNum(tmp) + 1;
    end
end
%% Sonar dataset
if strcmpi(select, 'sonar')
    %% load the sonar.txt
    disp("------------------------------------------------");
    disp("Using the Sonar dataset to test");

    data = importfile_sonar('sonar.txt', 1, 208);
    % simple transform: R(Rock) --> 1; M(Mine) --> 2;
    if ~isempty(data)
        [N, tt] = size(data);
    end 
    tmp_data = zeros(N, tt - 1); class_vec = zeros(N,1);
    for n = 1:N
        tclass = data{n,tt};
        if strcmpi(tclass, 'R')
            tmp_data(n,1) = 1;
            class_vec(n) = 1;
        elseif strcmpi(tclass, 'M')
            tmp_data(n,1) = 2;
            class_vec(n) = 2;
        end
        % other data
        for m = 1:tt-1
            tmp_data(n,m) = data{n,m};
        end        
    end
    
    data = [class_vec, tmp_data];
    if ~isempty(data) %^update
        [N, tt] = size(data);
    end
    
    Num_Attr  = tt-1; % The first is class;
    
    typeClass = [{'Rock'}; {'Mine'}];Num_class = length(typeClass);
    [ evid_set ] = loadSonar_new( data, typeClass );
    
    % re-order
%     s1 = data(data(:,1) == 1, :);
%     s2 = data(data(:,1) == 2, :);
%     data = [s1;s2];
    
    TotalNum = zeros(Num_class, 1);
    for n = 1:N
        tmp = data(n,1);
        TotalNum(tmp) = TotalNum(tmp) + 1;
    end
end
%% ionosphere dataset
if strcmpi(select, 'ionosphere')
    %% load the ionosphere.data
    disp("------------------------------------------------");
    disp("Using the ionosphere dataset to test");

    data = importfile_ionosphere('ionosphere.data', 1, 351);
    % simple transform: g(good) --> 1; b(bad) --> 2;
    if ~isempty(data)
        [N, tt] = size(data);
    end 
   tmp_data = zeros(N, tt - 1); class_vec = zeros(N,1);
    for n = 1:N
        tclass = data{n,tt};
        if strcmpi(tclass, 'g')
            tmp_data(n,1) = 1;
            class_vec(n) = 1;
        elseif strcmpi(tclass, 'b')
            tmp_data(n,1) = 2;
            class_vec(n) = 2;
        end
        % other data
        for m = 1:tt-1
            tmp_data(n,m) = data{n,m};
        end        
    end
    
    data = [class_vec, tmp_data];
        % re-order
    s1 = data(data(:,1) == 1, :);
    s2 = data(data(:,1) == 2, :);
    data = [s1;s2];
    
    if ~isempty(data) %^update
        [N, tt] = size(data);
    end   
    Num_Attr  = tt-1; % The first is class;
    
    typeClass = [{'good'}; {'bad'}];Num_class = length(typeClass);
    [ evid_set ] = loadWine_new( data, typeClass );
    TotalNum = zeros(Num_class, 1);
    for n = 1:N
        tmp = data(n,1);
        TotalNum(tmp) = TotalNum(tmp) + 1;
    end
end

%% Preparation of training -- k fold validation
format long;
if debug_kfoldValidation
    k = 3;
    foldSizes = computeFoldSizes(TotalNum, k);
    [X_sorted, y_sorted] = randSortAndGroup(data(:,2:Num_Attr+1), data(:,1), 1:Num_class);
    
    disp("------------------------------------------------");
    fprintf("[k = %d] k fold validation\n", k);
    accuracy_vec = zeros(k, 1);
    for roundNumber = 1 : k
        % Select the vectors to use for training and cross validation.
        X_sorted = [X_sorted,(1:N)'];        
        [X_train, y_train, X_val, y_val] = getFoldVectors(X_sorted, y_sorted, 1:Num_class, TotalNum, foldSizes, roundNumber);
        train_order = X_train(:,Num_Attr+1); 
        test_order = (1:N)';test_order(train_order) = [];
        % Train the classifier on the training set, X_train y_train
        trainset_all = evid_set(train_order);
        for n = 1:Num_class
            trainset{n} = trainset_all(y_train == n);   
%             tt = ones(numel(trainset{n}),1);
%             cred_vec{n} = tt/sum(tt);
            
            cred_vec{n} = credible_degree(trainset{n});
        end
        
        % Measure the classification accuracy on the validation set.
        Num_test = numel(test_order);tmp_type = zeros(Num_test, 1);
        resArray = cell(Num_test,1); realArray = cell(Num_test,1);
        maxScore = zeros(Num_test, 1);
        for n = 1:Num_test
            test_sample = evid_set{test_order(n)};
            cond_vec = zeros(Num_class, 1);
            for m = 1:Num_class
                InterDim_vec = interactionDimension(trainset{m}, test_sample);
                cond_vec(m) = sum(cred_vec{m}'.*InterDim_vec);
            end
            res = cond_vec == max(cond_vec);
            if sum(res)>1 % handle the situation of multiple maximum value(select the first)
                tmp = false(length(res), 1);
                for t = 1:length(res)
                    if res(t) == 1
                        tmp(t) = true;
                        break;
                    end
                end
                res = tmp;
            end
            tt = (1:Num_class);
            maxScore(n) = tt(res);
            %             resArray(n) = typeClass(res);realArray(n) = typeClass(data(test_order(n),1));
        end

        % report the result
        tt = abs(maxScore - y_val);
        Num_error = numel(find(tt>0));
        
        accuracy_vec(roundNumber) = (Num_test - Num_error)/Num_test;
        fprintf("[roundNumber = %d]:finished\n",roundNumber);pause(0.05);
    end
%     disp(sprintf("Accuracy : %f(%d/%d); Error: %f(%d/%d);",(Num_test - Num_error)/Num_test, (Num_test - Num_error),Num_test, Num_error/Num_test, Num_error, Num_test));
    mean_accuracy = mean(accuracy_vec); var_accuracy = std(accuracy_vec);
    fprintf("[Accuracy] mean-value: %f, variance: %f\n", mean_accuracy, var_accuracy);
    return;
end

%% using the new method to classify
N_trainset = 30; N_testset = TotalNum(1) - N_trainset;
train_order = 1:N_trainset;      %randperm(N_trainset);
test_order  = 1:(sum(TotalNum)); tmp_order = []; cumuValue = 0;
% get the train set 
for n = 1:Num_class
    trainset{n} = evid_set(train_order);  % randperm(N_trainset)+50;
    tmp_order = [tmp_order, train_order]; % for Class 1~2~3:{'setosa'}{'versicolor'}{'virginica'}
    
    cumuValue = cumuValue + TotalNum(n);
    train_order = (1:N_trainset)+cumuValue;
end
% train_order = (1:N_trainset)+50; %randperm(N_trainset)+50;
% trainset{2} = evid_set(train_order);% for Class 2:{'versicolor'}
% tmp_order = [tmp_order, train_order];
% 
% train_order = (1:N_trainset)+100;%randperm(N_trainset)+100;
% trainset{3} = evid_set(train_order);% for Class 3:
% tmp_order = [tmp_order, train_order];

test_order(tmp_order) = [];
% calculate the credible degree
cred_vec = zeros(Num_class, N_trainset);
for n = 1:Num_class
    cred_vec(n,:) = credible_degree(trainset{n});disp(sprintf('[Finished] training type %d', n));
end

Num_test = length(test_order);tmp_type = zeros(Num_test, 1);
resArray = cell(Num_test,1); realArray = cell(Num_test,1);
for n = 1:Num_test
    test_sample = evid_set{test_order(n)};
    disp('！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！！');
    
    cond_vec = zeros(Num_class, 1);
    for m = 1:Num_class
        InterDim_vec = interactionDimension(trainset{m}, test_sample);
        cond_vec(m) = sum(cred_vec(m,:).*InterDim_vec);
        % debug information
%         disp(strcat('[',num2str(test_order(n)),']:',sprintf('Finish the Class %d', m)));
    end
    
    
%     InterDim_vec = interactionDimension(trainset_c2, test_sample);
%     cond_vec_2 = sum(cred_vec_c2'.*InterDim_vec);
%     disp(strcat('[',num2str(test_order(n)),']:','Finish the Class 2'));
%     
%     InterDim_vec = interactionDimension(trainset_c3, test_sample);
%     cond_vec_3 = sum(cred_vec_c3'.*InterDim_vec);
%     disp(strcat('[',num2str(test_order(n)),']:','Finish the Class 3'));
    
    % judge the sample type
%     cond_vec = [cond_vec_1, cond_vec_2, cond_vec_3];
    res = cond_vec == max(cond_vec);
    if sum(res)>1 % handle the situation of multiple maximum value(select the first)
        tmp = false(length(res), 1);
        for t = 1:length(res)
            if res(t) == 1
                tmp(t) = true;
                break;
            end
        end
        res = tmp;
    end
    res_type = typeClass{res};
    resArray(n) = typeClass(res);realArray(n) = typeClass(data(test_order(n),1));
    
%     tmp_type(n) = find(cond_vec == max(cond_vec));
%     res_type_vec(n) = res_type;
    disp(strcat('[',num2str(test_order(n)),']:', typeClass{data(test_order(n),1)}, ' is belong to --->>',res_type));

end
%% Statistics
matchVec = zeros(Num_test,1);
for k = 1:Num_test % for each sample
    tmp = typeClass{data(test_order(k),1)};
    
    if strcmpi(resArray{k}, tmp)
            matchVec(k) = 1;
        else
            matchVec(k) = 0;
    end
end
disp("------------------------------------------------");
disp(sprintf("Accuracy : %f(%d/%d); Error: %f(%d/%d);",sum(matchVec)/Num_test, sum(matchVec),Num_test, (Num_test-sum(matchVec))/Num_test, Num_test-sum(matchVec), Num_test));