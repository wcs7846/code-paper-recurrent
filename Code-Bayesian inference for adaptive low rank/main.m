%% This script is designed to implement the method from the reference article.
% reference: Bayesian inference for adaptive low rank and sparse matrix estimation
% Author: MarkLHF £¨2751867750@qq.com£©
% Date: 2019/3/6
close all;clear;
detection = 0;
addpath(genpath('D:\Matlab_program\reference\TFOCS'));
%% IPI model
% input test image
ii = 1;
I = imread(strcat(num2str((ii)),'.bmp'));
figure;imshow(I);
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
I_double = im2double(I);
D = gen_patch_img(I_double, patchSize, slideStep);


[mm, nn]=size(D);
% parameters setting
% lambda =1/(sqrt(min(mm,nn)))*2;
% mu=sqrt(2*max(mm,nn))*4;     % this value is related to the degree of noise
% 
% SparseRatio = 0.05;
% E_true = zeros(size(D));
% Eomega = randsample(mm*nn, round(mm*nn*SparseRatio));
% E_true(Eomega) = rand(length(Eomega),1); % sparse component

%% Run ARLLR
tic
% [X_hat_t, ~, ~, E_hat_t] = VBRPCA(D,options);
% [X_hat, E_hat] = arllr(D, 10^-4);
[X_hat, E_hat, list_e_norm] = arllre(D, 100, 1);
toc
%% Reconstruct background and target image
[rstT, rstB] = res_patch_img_mean(E_hat, X_hat, I_double, patchSize, slideStep);
% Show the result
figure;
% Tips: the negative value is ignored
subplot(121),imshow(rstT .* (rstT>0), []),title('Target');
subplot(122),imshow(rstB .* (rstB>0), []),title('Background');
figure;mesh(rstT.* (rstT>0));title('residual');
% figure;plot(list_e_norm,'r.-');hold on;title('the norm of residual');
%% Run toolbox TFOCS
% Ref:http://cvxr.com/tfocs/
sigma = 1/(5*std(D(:)));
tao = 1.2*sqrt(2)*sigma^2;
lamda = 1;
tol = 5e-3;

% x_recon = tfocs(smooth_quad, {eye(mm), -D}, proj_l1(tao));

%% Detection part
if detection
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
    figure;imshow(I)%,'border','tight');%title('detect result');%,'border','tight'
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
end
% clear the addpath
rmpath(genpath('D:\Matlab_program\reference\TFOCS'));