%% This script is designed to implement the measure from the reference article.
% reference: Robust infrared small target detection using local steering kernel reconstruction
% Copyright:2018-9-4 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
clc;close all;clear;

% [Tips]: the path is data path, you need to CHANGE the path by yourself
addpath('D:\Matlab_program\self-code\test-data');
% step 1: prepare the test data
Img = imread('0002.bmp');
figure;imshow(Img);title('test image');
% step 2: get the patch image
s = 5; % must be odd number
[row, col, channel] = size(Img);
showImg = imresize(Img,[row,col]);
if channel > 1
    Img = rgb2gray(Img);
end
padding = (s-1)/2+s;
% 2.1: prepare the deta_coordinate matrix
[ coordinate, ~, len ] = pointVector(-(s-1)/2,(s-1)/2,-(s-1)/2,(s-1)/2);
% 2.2: padding
Img_ex = double(padarray(Img, [padding, padding], 'replicate'));
% 2.3: extract the patch image from the original image
LSK_patch = zeros(s,s,9);
lacm = zeros(row,col); % Local adaptive contrast measure based on regularized LSK
for nrow = 1:row
    for ncol = 1:col
        % patch image
        patch = Img_ex(nrow:nrow+3*s-1,ncol:ncol+3*s-1);
        % divide patch Image into 9 parts(3*3=9)
        % Tips: I calculate the gradient matrix before segmenting it
        [gx, gy] = derivative7(patch, 'x', 'y');
        gx = gx/(max(gx(:))+eps);
        gy = gy/(max(gy(:))+eps);
        for t = 1:9
            % Tips: The reason to set the tt and tp temporal parameter is at
            % the bottom
            tt = floor((t-1)/3);
            tp = mod(t-1,3);
            temp_gradient_x = gx(1+tt*s:(tt+1)*s, 1+tp*s:(tp+1)*s);
            temp_gradient_y = gy(1+tt*s:(tt+1)*s, 1+tp*s:(tp+1)*s);
            % Attention: I DONNOT UNDERSTAND WHY I NEED TO USE THE K 
%             win = (s-1)/2;
%             K = fspecial('disk', win);
%             K = K ./ K(win+1, win+1);
%             temp_gradient_x = temp_gradient_x .* K;
%             temp_gradient_y = temp_gradient_y .* K;
            % calculate the LSK of the patch
            % Tips: LSK must be greater than 0
            temp_lsk = LSK(coordinate, s, s, temp_gradient_x, temp_gradient_y);
            var = std(temp_lsk(:)); mean_v = mean(temp_lsk(:));
            if  var > eps
                temp_lsk_norm = normalise(temp_lsk, mean_v, var);
            else
                temp_lsk_norm = mean_v*ones(size(temp_lsk));
            end
            LSK_patch(:,:,t) = temp_lsk_norm;
        end
        % step 3: construct the fc and fN vector
        fc = reshape(LSK_patch(:,:,5), 1, []); % the central patch is 5th
        fN = zeros(8, length(fc)); tmp = 1;
        for t = 1:9
            if t ~= 5
                fN(tmp,:) = reshape(LSK_patch(:,:,t), 1, []);
                tmp = tmp + 1;
            end
        end
        % step 4: calculate the coefficient vector w of fN
        % --- model: Linear weighted model
        lamda = 1.0;
        w = pinv(fN*fN'+1.0*eye(8))*(fN*fc'); % pinv() == MP-inverse matrix(N*1)
        fT = w'*fN;
        lacm(nrow,ncol) = norm(fc-fT,2);
    end
end
figure;imshow(lacm,[]);title('LACM-LSK result');
figure;pcolor(lacm);title('LACM-LSK');colorbar;
figure;mesh(lacm);title('LACM-LSK');xlim([0,col]);ylim([0,row]);
%% Detection part
% binaryzation
mean_confidence = mean(lacm(:));
std_confidence  = std(lacm(:));
level = mean_confidence + 10*std_confidence;
bw = zeros(row,col);
bw(lacm>level) = 1;
bw = logical(bw);
% calculate the position
bw_img = bwlabel(bw);
stats = regionprops(bw_img);
Ar = cat(1, stats.Area);
target = stats(find(Ar ==max(Ar)));
% show 
pp = s;
figure;imshow(showImg);title('detect result');%,'border','tight'
for n = 1:length(target)
    pos = [(target(n).Centroid - pp/2) [pp pp]];
    rectangle('Position',1*pos,'EdgeColor','r');
end
%% the explaination for why set temporal parameter: tt and tp
% an analysis of 3 orders
% t  = 1,2,3,4,5,6,7,8,9
% tt = 0,0,0,1,1,1,2,2,2
% tp = 0,1,2,0,1,2,0,1,2
% So: tt = floor((t-1)/3); tp = mod((t-1),3);