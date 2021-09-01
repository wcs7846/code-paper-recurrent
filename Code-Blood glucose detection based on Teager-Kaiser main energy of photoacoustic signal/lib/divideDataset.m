function [ trainset, testset ] = divideDataset( dataset, proprety_list, size_train_set)
%DIVIDEDATASET
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  dataset        --> the dataset(struct)
         proprety_list  --> the list of proprety_list(dim:1D)
         size_train_set --> the size of train set
 Output: trainset --> the dataset of train (struct)
         testset  --> the dataset of test (struct)
 Tips:
 dataset.original: the original data of signal
 dataset.cf      : the center frequence proprety
 dataset.rmsf    : the root mean of frequence proprety
 dataset.ibw     : the instantaneous bandwidth proprety
 dataset.tk      : the TK energy proprety
%}
order = [1:2:length(proprety_list), 2:2:length(proprety_list)];
% order = randperm(length(proprety_list));
min_c = min(proprety_list); max_c = max(proprety_list);
norm_proprety_list = (proprety_list - min(proprety_list))/(max_c - min_c);
% norm_proprety_list = (proprety_list - mean(proprety_list))/(max_c - min_c)*2;

% size_train_set = 10;
size_test_set  = length(proprety_list) - size_train_set;
train_order = sort(order(1:size_train_set));
test_order  = sort(order(size_test_set+1:length(proprety_list)));
% train_order = 1:2:length(concentration_list);
% test_order  = 2:2:length(concentration_list);

if isfield(dataset, 'original')
    trainset.original = dataset.original(train_order, :);
    testset.original = dataset.original(test_order, :);
else
    disp('[warning]: dataset not have original section! ')
end
if isfield(dataset, 'cf')
    trainset.cf   = dataset.cf(train_order,:);
    testset.cf   = dataset.cf(test_order,:);
else
    disp('[warning]: dataset not have center frequence section! ')
end
if isfield(dataset, 'rmsf')
    trainset.rmsf = dataset.rmsf(train_order,:);
    testset.rmsf = dataset.rmsf(test_order,:);
else
    disp('[warning]: dataset not have root mean square frequence section! ')
end
if isfield(dataset, 'ibw')
    trainset.ibw  = dataset.ibw(train_order,:);
    testset.ibw  = dataset.ibw(test_order,:);
else
    disp('[warning]: dataset not have instantaneous bandwidth section! ')
end
if isfield(dataset, 'tk')
    trainset.tk  = dataset.tk(test_order,:);
    testset.tk  = dataset.tk(test_order,:);
else
    disp('[warning]: dataset not have tk energy section! ')
end
if isfield(dataset, 'fag')
    trainset.fag  = dataset.fag(test_order,:);
    testset.fag  = dataset.fag(test_order,:);
else
    disp('[warning]: dataset not have frequency attenuation gradient section! ')
end

trainset.concentraion = norm_proprety_list(train_order);% the list of concentration
testset.concentraion = norm_proprety_list(test_order);% the list of concentration
end

