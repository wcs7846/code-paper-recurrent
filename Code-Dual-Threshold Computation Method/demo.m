%% This script is designed to implement the method from the reference article.
% reference: Bowstring-Based Dual-Threshold Computation Method for Adaptive Canny Edge Detector
% Copyright:2019-3-1 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
clc;close all;clear;

Img = imread('4.bmp');
figure;imshow(Img);title('test image');
% Check: is color image, or not.
[row, col, dim] = size(Img);
if dim == 3
    Img = rgb2gray(Img);
end
Img = im2double(Img);
%% calculate the two threshold 
[ upper_thres, lower_thres ] = dualThreshold( Img, 1 );

%% show the result(but I don't know how to use two threshold)
fprintf("the high-threshold is %f\n", upper_thres);
fprintf("the low-threshold is %f\n",  lower_thres);

