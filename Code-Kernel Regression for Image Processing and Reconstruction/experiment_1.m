%% This script is designed to show the shape of LSK 
% and test this algorithm's performance when processing a big image
% Copyright:2018-10-28 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
clc;close all;clear;

% [Tips]: This script is just used as a test platform to comfire my ideas.
%         and I won't write enough comment to explain these operation
addpath('D:\Matlab_program\self-code\test-data');
addpath('.\support\')
Img = imread('0001.bmp');
figure;imshow(Img);title('test image');
% transform rgb to gray
[row, col, channel] = size(Img);
if channel > 1
    Img = rgb2gray(Img);
end
Img = im2double(Img);

% iteartive steering kernel regression (second order)
IT = 15;     % the total number of iterations
wsize = 11;  % the size of the local orientation analysis window
lambda = 1;  % the regularization for the elongation parameter
alpha = 0.5; % the structure sensitive parameter
h = 2.4;     % the global smoothing parameter
ksize = 5;   % the kernel size
r = 1;       % the upscaling factor
z = zeros(row, col, IT+1);
zx1 = zeros(row, col, IT+1);
zx2 = zeros(row, col, IT+1);
rmse = zeros(IT+1, 1);
z(:,:,1) = Img;
[zx1c, zx2c] = derivative7(Img, 'x', 'y');  
zx1(:,:,1) = zx1c;
zx2(:,:,1) = zx2c;
C = zeros(2, 2, row, col);

for i = 2 : IT+1
    % compute steering matrix
    C = steering(zx1(:,:,i-1), zx2(:,:,i-1), ones(size(Img)), wsize, lambda, alpha);
    % steering kernel regression
    [zs, zx1s, zx2s] = skr2_regular(z(:,:,i-1), h, C, r, ksize);
    z(:,:,i) = zs;
    zx1(:,:,i) = zx1s;
    zx2(:,:,i) = zx2s;  
    % display the No.iteration
    fprintf('[%d]Iteration\n',i);
    %figure(99); imagesc(zs); colormap(gray); axis image; pause(1);
end
% 
H11 = zeros(row, col);
H12 = zeros(row, col);
H22 = zeros(row, col);
for nrow = 1:row
    for ncol = 1:col
        H11(nrow, ncol) = C(1,1,nrow,ncol);
        H12(nrow, ncol) = C(1,2,nrow,ncol);
        H22(nrow, ncol) = C(2,2,nrow,ncol);
    end
end

% show
figure; imagesc(Img); colormap(gray); axis image;title('The noisy image, STD=25');
figure; imagesc(H11); colormap(gray); axis image;
figure; mesh(H11); title('H(1,1)');
figure; imagesc(H12); colormap(gray); axis image;
figure; mesh(H12); title('H(1,2)');