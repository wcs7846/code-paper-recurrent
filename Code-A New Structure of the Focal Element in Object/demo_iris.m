%% Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
% This script is designed to test the information fractals for evidential
% pattern classification
% Author: MarkLHF
% Date: 2020-1-7(first version)
%       2020-2-3
%% Experiment 1: Constructive database test
% Ref: Information Fractals for Evidential  pattern Classification
clc;clear all;close all;

addpath('./lib')
debug_showdistribution = 0;
%% use parpool to quicken the calculation
% parpool;

%% load the iris.csv
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
dataBase = []; TotalNum = zeros(Num_class, 1);
for n = 1:N
    tmp = struct('order', data{1}(n), 'Sepal_len', data{2}(n), 'Sepal_wid', data{3}(n), ...
        'Petal_len', data{4}(n), 'Petal_wid', data{5}(n), 'Species', data{6}(n));
    dataBase = [dataBase; tmp];
    %
    flag = false(Num_class, 1);
    for m = 1:Num_class
        if strcmpi(typeClass{m}, data{6}{n})
            TotalNum(m) = TotalNum(m) + 1;
        end
    end
end
%% create the full library(transform the source data to BPA)
% calculate the mean and std of every attributes of every class
mean_mat = zeros(Num_class, Num_Attr); % row: same class; col: same attribute
std_mat  = zeros(Num_class, Num_Attr);
record_mat = cell(Num_class, Num_Attr);
for m = 1:Num_Attr % for each attribute
    
    for k = 1:Num_class
        tmp_vec = zeros(TotalNum(k),1); t=1;
        for n = 1:N % for each sample
            switch(m)
                case 1
                    tmp_attr = dataBase(n).Sepal_len;
                case 2
                    tmp_attr = dataBase(n).Sepal_wid;
                case 3
                    tmp_attr = dataBase(n).Petal_len;
                case 4
                    tmp_attr = dataBase(n).Petal_wid;
            end
            if strcmpi(typeClass{k}, dataBase(n).Species)
                tmp_vec(t) = tmp_attr;
                t = t+1;
            end
        end
        record_mat(k,m) = {tmp_vec};
        % calcualte the mean and std value
        mean_mat(k,m) = mean(tmp_vec);
        std_mat(k,m) = std(tmp_vec);
    end
end

%% generator the BPA for each samples
evid_set = cell(N, 1);
for k = 1:N % for each sample
    tmp_samples = dataBase(k);
    tmp_vec = zeros(Num_Attr, 1);
    fe_vec = cell(Num_class, Num_Attr); % each row is a focal elements
    fe_bpa = zeros(Num_class, Num_Attr);
    for m = 1:Num_Attr
        % generate the BPA for each attributes
        switch(m)
            case 1
                tmp_attr = tmp_samples.Sepal_len;
            case 2
                tmp_attr = tmp_samples.Sepal_wid;
            case 3
                tmp_attr = tmp_samples.Petal_len;
            case 4
                tmp_attr = tmp_samples.Petal_wid;
        end
        tmp_bpa_vec = zeros(Num_class, 1);
        for n = 1:Num_class
            mu = mean_mat(n,m); s = std_mat(n,m);
            tmp_func = @(x)(1/sqrt(2*pi*s^2)*exp(-(x-mu)^2/(2*s^2)));
            
            tt = tmp_func(tmp_attr);
            if (tmp_attr<mu-3*s)||(tmp_attr>mu+3*s)
                tt = 0;
            end
            tmp_bpa_vec(n) = tt;
        end
        tmp_bpa_vec = tmp_bpa_vec/(sum(tmp_bpa_vec)+eps); % for No.44
        tt = [tmp_bpa_vec, (1:Num_class)'];
        
        % allocate the basic element
        tt2 = sortrows(tt, 1, 'descend');
        for n = 1:Num_class
            fe_vec(n,m) = typeClass(tt2(n,2));
            fe_bpa(n,m) = tt2(n,1);
        end
    end
    % construct the focal element
    tmp_evid = [];tt = 0;uni_set = createKnowledgeLib(typeClass, Num_Attr);
    for n = 1:Num_class
        tmp_bpa = mean(fe_bpa(n,:));
        tmp_fe.basicElement = {fe_vec(n,:)};
        tmp_fe.bpa = tmp_bpa;
        tmp_evid = [tmp_evid, tmp_fe]; tt = tt + tmp_bpa;
    end
    % for universal set
    tmp_fe.basicElement = uni_set;
    if tt > 1
        for n = 1:Num_class
            tmp_evid(n).bpa = tmp_evid(n).bpa/tt;
        end
        tmp_fe.bpa = 0;
    else
        tmp_fe.bpa = 1-tt;
    end
    tmp_evid = [tmp_evid, tmp_fe];
    evid_set(k) = {tmp_evid};
end
%% show the distribution of SL, SW, PL, PW
if debug_showdistribution
    for m = 1:Num_Attr % for each attribute
        figure;
        
        % draw the legend
        for n = 1:Num_class
            mu = mean_mat(n,m); s = std_mat(n,m);
            tmp_func = @(x)(1/sqrt(2*pi*s^2)*exp(-(x-mu)^2/(2*s^2)));
            switch(n)
                case 1
                    fplot(tmp_func, 'b--');hold on;
                case 2
                    fplot(tmp_func, 'r--');hold on;
                case 3
                    fplot(tmp_func, 'g--');hold on;
            end
        end
        % draw the scatter point
        for n = 1:Num_class
            record_tmp = record_mat{n,m};
            mu = mean_mat(n,m); s = std_mat(n,m);
            tmp_func = @(x)(1/sqrt(2*pi*s^2)*exp(-(x-mu).^2./(2*s^2)));
            y = tmp_func(record_tmp);
            switch(n)
                case 1
                    plot(record_tmp,y, 'b*');hold on;
                case 2
                    plot(record_tmp,y, 'r*');hold on;
                case 3
                    plot(record_tmp,y, 'g*');hold on;
            end
        end
        legend([typeClass;typeClass]);
        
        switch(m)
            case 1
                xlim([2,9]);ylabel('f_n_o_r_m_a_l(SL)');   % SL
            case 2
                xlim([1,5.5]);ylabel('f_n_o_r_m_a_l(SW)');  % SW
            case 3
                xlim([0,9]);ylabel('f_n_o_r_m_a_l(PL)');   % PL
            case 4
                xlim([0,4]);ylabel('f_n_o_r_m_a_l(PW)');    % PW
        end
        xlabel('length / cm');
        pause(0.5);
    end
end

%% Use the BPA to make the pattern recogition
N_trainset = 30; N_testset = TotalNum(1) - N_trainset;
train_order = 1:N_trainset;      %randperm(N_trainset);
test_order  = 1:(sum(TotalNum)); tmp_order = [];
% get the train set 
trainset_c1 = evid_set(train_order);% for Class 1:{'setosa'}
tmp_order = [tmp_order, train_order];

train_order = (1:N_trainset)+50; %randperm(N_trainset)+50;
trainset_c2 = evid_set(train_order);% for Class 2:{'versicolor'}
tmp_order = [tmp_order, train_order];

train_order = (1:N_trainset)+100;%randperm(N_trainset)+100;
trainset_c3 = evid_set(train_order);% for Class 3:{'virginica'}
tmp_order = [tmp_order, train_order];

test_order(tmp_order) = [];
% calculate the credible degree
% 论文中用来DM的例子为evi_set中的 23，33，44
tt_set = evid_set([23,33,44]);
credible_degree(tt_set);


cred_vec_c1 = credible_degree(trainset_c1);disp('[Finished] training type 1');
cred_vec_c2 = credible_degree(trainset_c2);disp('[Finished] training type 2');
cred_vec_c3 = credible_degree(trainset_c3);disp('[Finished] training type 3');

Num_test = length(test_order);tmp_type = zeros(Num_test, 1);
for n = 1:Num_test
    test_sample = evid_set{test_order(n)};
    disp('————————————————————————————————');
    
    InterDim_vec = interactionDimension(trainset_c1, test_sample);
    cond_vec_1 = sum(cred_vec_c1'.*InterDim_vec);
    disp(strcat('[',num2str(test_order(n)),']:','Finish the Class 1'));
    
    InterDim_vec = interactionDimension(trainset_c2, test_sample);
    cond_vec_2 = sum(cred_vec_c2'.*InterDim_vec);
    disp(strcat('[',num2str(test_order(n)),']:','Finish the Class 2'));
    
    InterDim_vec = interactionDimension(trainset_c3, test_sample);
    cond_vec_3 = sum(cred_vec_c3'.*InterDim_vec);
    disp(strcat('[',num2str(test_order(n)),']:','Finish the Class 3'));
    
    % judge the sample type
    cond_vec = [cond_vec_1, cond_vec_2, cond_vec_3];
    res = cond_vec == max(cond_vec);
    res_type = typeClass(res);
    
    tmp_type(n) = find(cond_vec == max(cond_vec));
    res_type_vec(n) = res_type;
    disp(strcat('[',num2str(test_order(n)),']:','is belong to --->>',res_type{1}));
end
