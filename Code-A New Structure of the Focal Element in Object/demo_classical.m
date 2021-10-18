%% Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
% This script is designed to test the information fractals for evidential
% pattern classification
% Author: MarkLHF
% Date: 2020-1-7(first version)
%       2020-10-20
%       2020-11-19(used to test the classical method)
%% Experiment 1: Constructive database test
% traditional focal element + traditional classification method
clc;clear all;close all;

addpath('./lib')
debug_showdistribution = 0;
debug_kfoldValidation = 1;
%% use parpool to quicken the calculation
% parpool;
% 'iris'; 'wine'; 'wdbc'; 'sonar'; 'ionosphere'
select = 'ionosphere';
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
    
    data = tmp_data;
    if ~isempty(data) %^update
        [N, tt] = size(data);
    end
    
    Num_Attr  = tt-1; % The first is class;
    
    typeClass = [{'malignant'}; {'benign'}];
    Num_class = length(typeClass);
    [ evid_set ] = loadWDBC_new( data, typeClass );
    
    % re-order
    s1 = data(data(:,1) == 1, :);
    s2 = data(data(:,1) == 2, :);
    data = [s1;s2];
    
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
    if ~isempty(data) %^update
        [N, tt] = size(data);
    end
    
    Num_Attr  = tt-1; % The first is class;
    
    typeClass = [{'good'}; {'bad'}];Num_class = length(typeClass);
    [ evid_set ] = loadWine_new( data, typeClass );
    
    % re-order
    s1 = data(data(:,1) == 1, :);
    s2 = data(data(:,1) == 2, :);
    data = [s1;s2];
    
    TotalNum = zeros(Num_class, 1);
    for n = 1:N
        tmp = data(n,1);
        TotalNum(tmp) = TotalNum(tmp) + 1;
    end
end

%% Preparation of training -- k fold validation
if debug_kfoldValidation
    k = 3;
    foldSizes = computeFoldSizes(TotalNum, k);
    [X_sorted, y_sorted] = randSortAndGroup(data(:,2:Num_Attr+1), data(:,1), 1:Num_class);
    
    disp("------------------------------------------------");
    disp(sprintf("[k = %d] k fold validation", k));
    accuracy_vec = zeros(k, 1);
    for roundNumber = 1 : k
        % Select the vectors to use for training and cross validation.
        [X_train, y_train, X_val, y_val] = getFoldVectors(X_sorted, y_sorted, 1:Num_class, TotalNum, foldSizes, roundNumber);
        
        % Train the classifier on the training set, X_train y_train
        % 'svm' 'REPtree' 'knn'
        select_method = 'svm';
        if strcmpi(select_method, 'svm')
            % multiple-svm == M classifier.
            model_vec = cell(Num_class,1);
            [N_train,~] = size(X_train);
            for n = 1:numel(typeClass)
                group = zeros(N_train, 1);
                group(y_train == n) = 1;
                
                tmp_model = fitcsvm(X_train, group,'Standardize',true,'KernelFunction','RBF','KernelScale',2);
                
                model_vec(n) = {tmp_model};
            end
        elseif strcmpi(select_method, 'REPtree')
            tmp_model = fitctree(X_train, y_train);
        elseif strcmpi(select_method, 'KNN')
            tmp_model = fitcknn(X_train, y_train);
        end
        
        % Measure the classification accuracy on the validation set.
        Num_test = numel(y_val);
        if strcmpi(select_method, 'svm')
            Scores = zeros(numel(y_val),numel(typeClass));
            for n = 1:numel(typeClass)
                [~,score_vec] = predict(model_vec{n}, X_val);
                Scores(:,n) = score_vec(:,2); % using the second column
            end
            [~,maxScore] = max(Scores,[],2);
            % verity
        elseif strcmpi(select_method, 'REPtree')
            [maxScore] = predict(tmp_model, X_val);
        elseif strcmpi(select_method, 'knn')
            [maxScore] = predict(tmp_model, X_val);
        end
        % report the result
        tt = abs(maxScore - y_val);
        Num_error = numel(find(tt>0));
        
        accuracy_vec(roundNumber) = (Num_test - Num_error)/Num_test;
    end
%     disp(sprintf("Accuracy : %f(%d/%d); Error: %f(%d/%d);",(Num_test - Num_error)/Num_test, (Num_test - Num_error),Num_test, Num_error/Num_test, Num_error, Num_test));
    mean_accuracy = mean(accuracy_vec); var_accuracy = std(accuracy_vec);
    disp(sprintf("[Accuracy] mean-value: %f, variance: %f", mean_accuracy, var_accuracy));
    return;
end


%% Preparation of training --
N_trainset = 20; N_testset = TotalNum(1) - N_trainset;
train_order = 1:N_trainset;      %randperm(N_trainset);
test_order  = 1:(sum(TotalNum)); tmp_order = []; cumuValue = 0;
% get the train set
for n = 1:Num_class
    trainset{n} = data(train_order,:);  % randperm(N_trainset)+50;
    tmp_order = [tmp_order, train_order]; % for Class 1~2~3:{'setosa'}{'versicolor'}{'virginica'}
    
    cumuValue = cumuValue + TotalNum(n);
    train_order = (1:N_trainset)+cumuValue;
end
test_order(tmp_order) = [];
testset = data(test_order,:);
%% Train classical method
% 'svm' 'REPtree' 'knn'
select_method = 'svm';
if strcmpi(select_method, 'svm')
    % multiple-svm == M classifier.
    model_vec = cell(Num_class,1);
    for n = 1:numel(typeClass)
        tmp_data_is  = trainset{n};
        tmp_order = 1:Num_class; tmp_order(tmp_order == n) = [];tmp_data_not = [];
        for m = 1:length(tmp_order)
            tmp_data_not = [tmp_data_not; trainset{tmp_order(m)}];
        end
        ttmp_data = [tmp_data_is(:,2:(Num_Attr+1)); tmp_data_not(:, 2:(Num_Attr+1))];
        [len,~] = size(tmp_data_not);
        group = [tmp_data_is(:,1); zeros(len,1)];
        
        tmp_model = fitcsvm(ttmp_data, group,'Standardize',true,'KernelFunction','RBF','KernelScale',1);
        
        model_vec(n) = {tmp_model};
    end
elseif strcmpi(select_method, 'REPtree')
    ttmp_data = [];
    for n = 1:numel(trainset)
        ttmp_data = [ttmp_data; trainset{n}];
    end
    tmp_model = fitctree(ttmp_data(:,2:(Num_Attr+1)), ttmp_data(:,1));
elseif strcmpi(select_method, 'KNN')
    ttmp_data = [];
    for n = 1:numel(trainset)
        ttmp_data = [ttmp_data; trainset{n}];
    end
    tmp_model = fitcknn(ttmp_data(:,2:(Num_Attr+1)), ttmp_data(:,1));
end


%% predict
disp("------------------------------------------------");
Num_test = numel(test_order);
if strcmpi(select_method, 'svm')
    Scores = zeros(numel(test_order),numel(typeClass));
    for n = 1:numel(typeClass)
        [~,score_vec] = predict(model_vec{n},testset(:,2:Num_Attr+1));
        Scores(:,n) = score_vec(:,2); % using the second column
    end
    [~,maxScore] = max(Scores,[],2);
    % verity
elseif strcmpi(select_method, 'REPtree')
    [maxScore] = predict(tmp_model, testset(:,2:Num_Attr+1));
elseif strcmpi(select_method, 'knn')
    [maxScore] = predict(tmp_model, testset(:,2:Num_Attr+1));
end
% report the result
tt = abs(maxScore - testset(:,1));
Num_error = numel(find(tt>0));

disp(sprintf("Accuracy : %f(%d/%d); Error: %f(%d/%d);",(Num_test - Num_error)/Num_test, (Num_test - Num_error),Num_test, Num_error/Num_test, Num_error, Num_test));
