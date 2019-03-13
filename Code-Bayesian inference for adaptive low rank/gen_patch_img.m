function patchImg = gen_patch_img(img, patchSize, slideStep)
%% GEN_PATCH_IMG
% This matlab code generates the patch-image for infrared 
% patch-image model.
% Copyright:2019-3-12 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
%                     Tianfang Zhang,  UESTC.(e-mail:sparkcarleton@gmail.com)
%{
 detail
 Input:  img       --> the input image(2D, and uint8 type);
         patchSize --> the size of patch(size = patchSize*patchSize)
         slideStep --> the step of slide window
 Output: patchImg  --> the patch image(2D, double)
%}
if ~exist('patchSize', 'var')
    patchSize = 50;
end
if ~exist('slideStep', 'var')
    slideStep = 10;
end

% img = reshape(1:9, [3 3])
% img = reshape(1:12, [3 4])

[row, col] = size(img);
D = zeros(patchSize*patchSize, (length(1:slideStep:row-patchSize)*length(1:slideStep:col-patchSize)));
counter = 1;

%% arrayfun version, identical to the following for-loop version
for i = 1:slideStep:row-patchSize
    for j = 1:slideStep:col-patchSize
        tmp_patch = img(i:i+patchSize-1, j:j+patchSize-1);
        D(:, counter) = reshape(tmp_patch, patchSize*patchSize, 1);
        counter = counter + 1;
    end
end
%% output
patchImg = D;
%% for-loop version
% patchImg = zeros(patchSize * patchSize, rowPatchNum * colPatchNum);
% k = 0;
% for col = colPosArr
%     for row = rowPosArr
%         k = k + 1;
%         tmp_patch = img(row : row + patchSize - 1, col : col + patchSize - 1);
%         patchImg(:, k) = tmp_patch(:);
%     end
% end

