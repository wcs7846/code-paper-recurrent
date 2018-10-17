% This script is designed to implement the measure from the reference
% article.
% reference: An infrared small target detection method based on multiscale local homogeneity measure
% operation
% Copyright:2018-9-3 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
clc;close all;clear;

% [Tips]: the path is data path, you need to CHANGE the path by yourself
addpath('D:\Matlab_program\self-code\test-data');
% step 1: prepare the test data
Img = imread('2.bmp');
figure;imshow(Img);title('test image');

[row, col] = size(Img);
matrix_w    = zeros(row, col);
matrix_w_ex = zeros(row, col, 6-3+1);
% step 2: set multiscale 
t = 1;
for L = 3:2:9
    % 2.1£º add a padding to Img
    padding = (L-1)/2+L;
    Img_ex = impad(Img, padding, 'replicate');
    
    % step 3: use LHM to calculate the feature map of local homogeneity measure
    [matrix_w] = LHM(Img_ex, L);
    matrix_w_ex(:,:,t) = matrix_w;
    t = t+1;
end
% step 3: merge multiscale LHM to a 2D-matrix
result = zeros(row, col);
for nrow = 1:row
    for ncol = 1:col
        % calculate the max w
        Wvector = matrix_w_ex(nrow,ncol,:);
        result(nrow,ncol) = max(Wvector);
    end
end
figure;imshow(result,[]);title('LHM-2D');
figure;mesh(result);title('LHM-3D');
%% Detection part
% binaryzation
mean_confidence = mean(result(:));
std_confidence  = std(result(:));
level = mean_confidence + 10*std_confidence;
bw = zeros(row,col);
bw(result>level) = 1;
bw = logical(bw);
% calculate the position
bw_img = bwlabel(bw);
stats = regionprops(bw_img);
Ar = cat(1, stats.Area);
target = stats(find(Ar ==max(Ar)));
% show 
pp = 5;
pos = [(target.Centroid - pp/2) [pp pp]];
figure;imshow(imresize(Img,2*[row,col]),'border','tight');title('detect result');rectangle('Position',2*pos,'EdgeColor','r');
