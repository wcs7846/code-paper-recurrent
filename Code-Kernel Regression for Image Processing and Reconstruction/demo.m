%% This script is designed to implement the measure from the reference article.
% reference: Kernel Regression for Image Processing and Reconstruction
% Copyright:2018-10-26 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
clc;close all;clear;

% [Tips]: the path is data path, you need to CHANGE the path by yourself
addpath('D:\Matlab_program\self-code\test-data');
% step 1: prepare the test data
Img = imread('0001.bmp');
figure;imshow(Img);title('test image');
