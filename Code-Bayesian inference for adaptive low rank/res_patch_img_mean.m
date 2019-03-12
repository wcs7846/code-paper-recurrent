function [recImg_T, recImg_B] = res_patch_img_mean(patchImg_T, patchImg_B, img, patchSize, slideStep)
%% RES_PATCH_IMG_MEAN
% This matlab code recovers the image from patch-image using 
% mean filter for the sake of speed.
% Copyright:2019-3-12 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
%                     Tianfang Zhang,  UESTC.(e-mail:sparkcarleton@gmail.com)
%{
 detail
 Input:  patchImg_T --> the patch sparse   image(2D, and double type);
         patchImg_B --> the patch low-rank image(2D, and double type);
         img        --> the input image         (2D, and uint8 type);
         patchSize  --> the size of patch(size = patchSize*patchSize)
         slideStep  --> the step of slide window
 Output: recImg_T   --> the reconstruction sparse image(2D, double)
         recImg_B   --> the reconstruction low-rank image(2D, double)
%}
%% get the basic information and allocate the store space
[row, col] = size(img);

reconsT = zeros(row, col, 100);
reconsB = zeros(row, col, 100);
countMatrix = zeros(row, col);
index = 1;
%% reconstruct the original image
for i = 1:slideStep:row-patchSize
    for j = 1:slideStep:col-patchSize
        % Count the time each pixel used and record the value
        repatch_T = reshape(patchImg_T(:,index), patchSize, patchSize);
        repatch_B = reshape(patchImg_B(:,index), patchSize, patchSize);
        countMatrix(i:i+patchSize-1, j:j+patchSize-1) = countMatrix(i:i+patchSize-1, j:j+patchSize-1) + 1;
        
        % Record the value of each pixel
        for ii = i:i+patchSize-1
            for jj = j:j+patchSize-1
                reconsT(ii, jj, countMatrix(ii,jj)) = repatch_T(ii-i+1, jj-j+1);
                reconsB(ii, jj, countMatrix(ii,jj)) = repatch_B(ii-i+1, jj-j+1);
            end
        end
        
        index = index + 1;
    end
end
%% deal with the overlapping region
rstT = zeros(row, col);
rstB = double(img);

for i = 1:row
    for j = 1:col
        % Median
        if countMatrix(i ,j) > 0
            vectorT = reconsT(i, j, 1:countMatrix(i,j));
            vectorB = reconsB(i, j, 1:countMatrix(i,j));
            
            rstT(i, j) = mean(vectorT);
            rstB(i, j) = mean(vectorB);
        end
    end
end
%% output 
recImg_T = rstT;
recImg_B = rstB;