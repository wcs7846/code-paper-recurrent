function [ X, E ] = arllr( input_args, sita )
%% ARLLR 
% reference: Unsupervised ridge detection using second order anisotropic Gaussian kernels
% Copyright:2019-2-27 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
%{
 detail
 Input:  input_args --> the input image(2D, and double type,range[0,1]);
         sita       --> a parameter(10*-4)
 Output: X --> the reconstruction(2D, double)
         E --> the residual
%}
%% get the basic information and allocate the store space
% sita = 10^-4; 
[row, col] = size(input_args);
X = zeros(row, col);
E = zeros(row, col);

sigma = 1/(5*std(input_args(:)));
tao = 1.2*sqrt(2)*sigma^2;
total_itera = 1;

% get the singular value matrix of Y and X

[U, t_sig_Y, V] = svd(input_args, 'econ'); 
sig_Y = diag(t_sig_Y);

%% main loop
for k = 1:total_itera
    %% calculate the every singular
    sig_X = proximal_log(sig_Y, tao, sita, 1);
    % reconstruct
    S_x = diag(sig_X);  
    X = U * S_x * V';
    %% update
    E = input_args - X;
end


end

