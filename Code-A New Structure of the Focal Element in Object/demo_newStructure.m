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
end

%% traditional classification method
% For each sample, because there is only one evidence, so we make the
% decision by using PPT transform directly
resArray = cell(N,1); pArray = zeros(N,1);
for k = 1:N % for each sample
    % PPT transform
    tmp_evid = evid_set{k};
    
    pVec = zeros(1,Num_class);
    for m = 1:Num_class
        % the probability for each class
        N_fe = length(tmp_evid);
        for n = 1:N_fe
            tmp_fe = tmp_evid(n);
            tmp_bpa = tmp_evid(n).bpa;
            if isempty(tmp_fe.basicElement) % handle the [] set
                continue;
            end
            N_be = length(tmp_fe.basicElement{1}); % judge
            % compare
            tmp_vec = zeros(1,N_be);
            for l = 1:N_be
                tmp_be = tmp_fe.basicElement{1}(l);
                if strcmpi(tmp_be{1}, typeClass{m})
                    tmp_vec(l) = 1;
                else
                    tmp_vec(l) = 0;
                end
            end
            % calculate the prop
            pVec(m) = pVec(m) + tmp_bpa*(sum(tmp_vec)/N_be);
        end
    end
    % make a decision
    tN = length(pVec); maxP = pVec(1); res = typeClass(1);
    for n = 1:tN
        if pVec(n) > maxP
            maxP = pVec(n);
            res = typeClass(n);
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
disp(sprintf("Accuracy : %f(%d/%d); Error: %f(%d/%d);",sum(matchVec)/N, sum(matchVec),N, (N-sum(matchVec))/N, N-sum(matchVec), N));

