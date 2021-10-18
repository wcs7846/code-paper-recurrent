%% Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
% This script is designed to test the information fractals for evidential
% pattern classification
% Author: MarkLHF
% Date: 2020-1-7(first version)
%       2020-10-20
%% Experiment 1: Constructive database test
% traditional focal element + traditional classification method
clc;clear all;close all;

addpath('./lib')
debug_showdistribution = 0;
%% use parpool to quicken the calculation
% parpool;
% 'iris'; 'wine'; 'wdbc'; 'sonar'; 'ionosphere'
select = 'iris';
% Tips: the first data stands for the classification.
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
    evid_set = loadIRIS_tradition(data, TotalNum, typeClass);
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
    [ evid_set ] = loadWine_tradition( data, typeClass );
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
    [ evid_set ] = loadWDBC_tradition( data, typeClass );
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
    
    typeClass = [{'Rock'}; {'Mine'}];
    [ evid_set ] = loadSonar_tradition( data, typeClass );
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
    
    typeClass = [{'good'}; {'bad'}];
    [ evid_set ] = loadWine_tradition( data, typeClass );
end
%% evidence combination
resArray = cell(N,1); pArray = zeros(N,1); 
for k = 1:N % for each sample
    tmp_samples = evid_set{k};
    tmp_evid = tmp_samples{1};
    for m = 2:Num_Attr
        tmp_evid = Murphy_evidMerge_trad(tmp_evid, tmp_samples{m});
    end
%     k_mat(k,:) = k_vec;
    % pignistic function
    discernmentFrame.type = typeClass;
    p = pignistic_prop_trans(tmp_evid, discernmentFrame);
    
    % make a decision
    tN = length(p); maxP = p(1).p; res = p(1).type;
    for n = 1:tN
        if p(n).p > maxP
            maxP = p(n).p;
            res = p(n).type;
        end
    end
    % show and record the result
    resArray(k) = res;
    pArray(k) = maxP;
   
    disp(sprintf("[%d] our: %s, real: %s",k, res{1}, typeClass{data(k,1)}));pause(0.01);
end
%% Statistics
matchVec = zeros(N,1);
for k = 1:N % for each sample
    tmp = typeClass{data(k,1)};
    
    if strcmpi(resArray{k}, tmp)
            matchVec(k) = 1;
        else
            matchVec(k) = 0;
    end

end
disp("------------------------------------------------");
disp(sprintf("[dataset] Accuracy : %f(%d/%d); Error: %f(%d/%d);",sum(matchVec)/N, sum(matchVec),N, (N-sum(matchVec))/N, N-sum(matchVec), N));

% %% anlaysis
% true_K = k_mat(:,matchVec == 1);
% false_K = k_mat(:,abs(matchVec-1) == 1);
% 
% mean_true_k = mean(true_K, 1);
% mean_false_k = mean(false_K, 1);