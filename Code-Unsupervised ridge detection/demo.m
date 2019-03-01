%% This script is designed to implement the method from the reference article.
% reference: Unsupervised ridge detection using second order anisotropic Gaussian kernels
% Copyright:2019-2-27 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
clc;close all;clear;

% [Tips]: the path is data path, you need to CHANGE the path by yourself
% addpath('D:\Matlab_program\self-code\test-data');
% step 1: prepare the test data
% address = 'D:\Matlab_program\self-code\code-paper-recurrent\Code-Unsupervised ridge detection\Ghent_University_Fungal_Images_I\Real_Images\';
% Img = imread(strcat(address, 'Crop001.png')); % Crop001.png;
Img = imread('7.bmp');
figure;imshow(Img);title('test image');
figure;imshow(Img,[]);title('test image(enhance)');
% Check: is color image, or not.
[row, col, dim] = size(Img);
if dim == 3
    Img = rgb2gray(Img);
end
%% <1> generate kernel created from some parameter(D, S, A)
% parameter instruction
% direction : D = (d1, d2, ... , dk) and di in [0 pi]
% scale     : S = (s1, s2, ... , sm) and si > 0
% anisotropy: A = (a1, a2, ... , an) and ai > 1
d_step = 20/180*pi;
D = 0:d_step:pi;
S = [1,2,4,6];
A = 1:0.1:1.5;
% generate kernel parameter matrix(3*N); N = the number of kernel
para_mat = getKernelParameter(D,S,A);
[~, N] = size(para_mat);

window_size = 11; % must be even
kernel_para = getAnisotropyKernel(para_mat, window_size);
% filtering
% second order derivative of an image by using AGK(anisotropic Gaussian kernels)
padding = window_size-1;
response_filter = zeros(row + padding, col + padding, N);
row_padd = row + padding; col_padd = col + padding;
for n = 1:N
    % extract the kernel
    temp_kernel = reshape(kernel_para(:,n), window_size, window_size);
    temp_para = para_mat(:,n);
    % d = direction; s = scale; a = anisotropy;
    di = temp_para(1);  si = temp_para(2);  ai = temp_para(3);
    [wx,wy] = meshgrid(1:window_size, 1:window_size);
    
    temp_kernel_coeff = ((ai/si)^2)*((wx.*cos(di)+wy.*(sin(di)^2))./(ai^-2*si^2)-1);
    temp_kernel = temp_kernel_coeff.*temp_kernel;
    % filter
    temp_response = conv2(Img, temp_kernel, 'full');

    % store the response
    response_filter(:,:,n) = temp_response;
%     figure;imshow(temp_response,[]);title('temp_response');
end 
%% select the maximum response for each pixel and save
% calculate the maximum response (propose: calculate the maximum intensity)
maxRes = max(response_filter, [], 3);
% locate the maximum response    (propose: calculate the optimal direction)
optDir = zeros(row_padd, col_padd);
optScale = zeros(row_padd, col_padd);
optAnisotropy = zeros(row_padd, col_padd);
for nrow = 1:row_padd
     for ncol = 1:col_padd
         maxValue = maxRes(nrow, ncol);
         temp_loc = find(response_filter(nrow, ncol, :) == maxValue);
         % check the uniqueness of maximum
         if length(temp_loc) > 1
             % [TODO]: calculate the optimal direction when the maximum
             % isnot unique
             temp_loc = temp_loc(1);
         end
         % calculate the optimal direction
         optimal_para = para_mat(:,temp_loc);
         
         optDir(nrow, ncol) = optimal_para(1);
         optScale(nrow, ncol) = optimal_para(2);
         optAnisotropy(nrow, ncol) = optimal_para(3);
     end
end
%% determining which pixels have a local maximum intensity
% [Tips]£º this part refer to David Young's code in MathWorks
% calculate the gradients
[G, REGION] = gradients_n(maxRes);
% non-maximum suppression
[Response_suppress, gMag] = nonmaxSuppress(G);
% ensure the same size between the non-maximum suppression process before
% and after
Response_suppress = padarray(Response_suppress, [1, 1], 'replicate');
gMag = padarray(gMag, [1, 1], 'replicate'); % gMag: the magnitude of gradient
%% determing threshold using OTSU(the part is different than this paper)
% but the result is same
otsu_level = graythresh(gMag);

%% applying hysteresis to generate the binary, thin map
result_detect_t = hystThresh(Response_suppress, gMag, otsu_level);
result_detect = result_detect_t(1:row, 1:col);
%% show the result
figure;imshow(maxRes(1:row, 1:col),[]);title('maxRes');
figure;imshow(Response_suppress(1:row, 1:col));title('Response suppress');
figure;imshow(result_detect);    title('result detect');

%% pure otsu
Img_double = im2double(Img);
otsu_level = graythresh(Img_double);
otsu_segment = zeros(size(Img));
otsu_segment(Img_double > otsu_level) = 1;
figure;imshow(otsu_segment); title('OTSU segment');