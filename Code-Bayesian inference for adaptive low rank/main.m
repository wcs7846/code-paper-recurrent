%% This script is designed to implement the method from the reference article.
% reference: Bayesian inference for adaptive low rank and sparse matrix estimation
% Author: MarkLHF £¨2751867750@qq.com£©
% Date: 2019/3/6
close all;clear;

%% IPI model
% input test image
ii = 30;
I = imread(strcat(num2str((ii)),'.bmp'));
[p, q, z]=size(I);

% rgb2gray
if z==3
    I=rgb2gray(I);
end

% set patch size (adjustable)
patchSize = 50;
slideStep = 10;

% patch model
% I = im2double(I);
D = gen_patch_img(I,patchSize, slideStep);
D = im2double(D);

[mm, nn]=size(D);
% parameters setting
% lambda =1/(sqrt(min(mm,nn)))*2;
% mu=sqrt(2*max(mm,nn))*4;     % this value is related to the degree of noise

SparseRatio = 0.05;
E_true = zeros(size(D));
% Eomega = randsample(mm*nn, round(mm*nn*SparseRatio));
% E_true(Eomega) = rand(length(Eomega),1); % sparse component

%% Setting VBRPCA parameter

%% Run VBRPCA
tic
% [X_hat, A_hat, B_hat, E_hat] = VBRPCA(D,options);

toc

%% reverse background and target image
corrupted  = res_patch_img_mean(E_hat, I, patchSize, slideStep);  % target image
background = res_patch_img_mean(X_hat, I, patchSize, slideStep);
% corrupted = E_hat;
% background = X_hat;

figure;imshow(corrupted,[]);
figure;imshow(background,[]);
%% Detection part
% binaryzation
mean_confidence = mean(corrupted(:));
std_confidence  = std(corrupted(:));
level = mean_confidence + 5*std_confidence;
bw = zeros(p,q);
bw(corrupted>level) = 1;
bw = logical(bw);
% calculate the position
bw_img = bwlabel(bw);
stats = regionprops(bw_img, 'Area'); %#ok<MRPBW>
Ar = cat(1, stats.Area);
target = stats(find(Ar ==max(Ar))); %#ok<FNDSB>
% show 
pp = 10;
figure;imshow(I,'border','tight');%title('detect result');%,'border','tight'
for n = 1:length(target)
    pos = [(target(n).Centroid - pp/2) [pp pp]];
    rectangle('Position',1*pos,'EdgeColor','r');
end
%% save the result to BMP image
% save target result
save_result = corrupted - min(corrupted(:));
save_result = save_result/max(save_result(:));
save_result_tmp = uint8(save_result*255);
% imwrite(save_result_tmp, strcat('target_',num2str((ii)),'.bmp'));
% save background result
save_result = background - min(background(:));
save_result = save_result/max(save_result(:));
save_result_tmp = uint8(save_result*255);
% imwrite(save_result_tmp, strcat('background_',num2str((ii)),'.bmp'));
