function [ output_args ] = getAnisotropyKernel( para_mat, window_size )
%% GETANISOTROPYKERNEL 
% reference: Unsupervised ridge detection using second order anisotropic Gaussian kernels
% Copyright:2019-2-27 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
%{
 detail
 Input:  para_mat    --> the parameter matrix(3*N); N = the number of kernel
         window_size --> the size of window(ws) 
 Output: output_args --> the kernel coefficient set((ws*ws)*N); N = the number of kernel
%}
%% get the basic information
[~, N] = size(para_mat);
[wx,wy] = meshgrid(1:window_size, 1:window_size);
output_args = zeros(window_size*window_size, N);
%% generate the kernel 
for n = 1:N
    % extract parameter
    % d = direction; s = scale; a = anisotropy;
    temp_para = para_mat(:,n);
    di = temp_para(1);  si = temp_para(2);  ai = temp_para(3);
    % calcualte rotation matrix
    rotMat_p = getRotationMatrix(di);  % positive
    rotMat_n = getRotationMatrix(-di); % negative
    aiMat = [ai.^2 0;0 ai.^2];
    
    wx_row = reshape(wx, 1, window_size*window_size);
    wy_row = reshape(wy, 1, window_size*window_size);
    
    wxy = [wx_row; wy_row];
    
    faiMat = diag(wxy' * rotMat_n * aiMat * rotMat_p * wxy);
    fai = reshape(faiMat, window_size, window_size);
    % calculate the kernel(size: window_size*window_size)
    temp_kernel = 1/(2*pi*si.^2)*exp(-fai./(2*si.^2));
    % store the kernel
    output_args(:, n) = reshape(temp_kernel, window_size*window_size, 1);
end
end

function  rotMat = getRotationMatrix(sita)
% this function is designed to generate the 2D rotation matrix
% [Tips]: sita's unit is radian(not degreee)
t_cos = cos(sita);
t_sin = sin(sita);

rotMat = [ t_cos t_sin;
          -t_sin t_cos;];
end

