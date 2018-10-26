%% This script is designed to implement the measure from the reference article.
% reference: Robust infrared small target detection using local steering kernel reconstruction
% Copyright:2018-9-4 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
clc;close all;clear;

% [Tips]: This script is just used as a test platform to comfire my ideas.
%         and I won't write enough comment to explain these operation
addpath('D:\Matlab_program\self-code\test-data');
Img = imread('0001.bmp');
figure;imshow(Img);title('test image');
% transform rgb to gray
[row, col, channel] = size(Img);
if channel > 1
    Img = rgb2gray(Img);
end
Img = im2double(Img);
% prepare the deta_coordinate matrix
[ coordinate, ~, len ] = pointVector(-(col-1)/2,(col-1)/2,-(row-1)/2,(row-1)/2);
% calculate the gradient matrix
[gx, gy] = derivative7(Img, 'x', 'y');
gx = gx/(max(gx(:))+eps);
gy = gy/(max(gy(:))+eps);
% calculate the LSK
LSK_temp = LSK(coordinate, row, col, gx, gy);
var = std(LSK_temp(:)); mean_v = mean(LSK_temp(:));
if  var > eps
    temp_lsk_norm = normalise(LSK_temp, mean_v, var);
else
    temp_lsk_norm = mean_v*ones(size(LSK_temp));
end
% show
figure;imshow(temp_lsk_norm,[]);