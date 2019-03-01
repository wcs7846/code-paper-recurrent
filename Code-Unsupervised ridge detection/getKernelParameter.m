function [ output_args ] = getKernelParameter( D, S, A)
%% GETKERNELPARAMETER 
% reference: Unsupervised ridge detection using second order anisotropic Gaussian kernels
% Copyright:2019-2-27 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
%{
 detail
 Input:  D = (d1, d2, ... , dk) and di in [0 pi](direction)  (1-D)
         S = (s1, s2, ... , sm) and si > 0      (scale)      (1-D)
         A = (a1, a2, ... , an) and ai > 1      (anisotropy) (1-D)
 Output: output_args --> the parameter matrix(3*N); N = the number of kernel
%}
%% allocate space to store the parameter
len_D = length(D);
len_S = length(S);
len_A = length(A);

N = len_A * len_D * len_S;

output_args = zeros(3, N);
%% generate the parameter matrix
% output_args(1,:) = D
% output_args(2,:) = S
% output_args(3,:) = A
% <1> set the D parameter
temp = repmat(D, len_S*len_A, 1);
output_args(1, :) = reshape(temp, [1, N]);
% <2> set the S parameter
temp = repmat(S, len_D*len_A, 1);
output_args(2, :) = reshape(temp, [1, N]);
% <3> set the A parameter
temp = repmat(A, len_S*len_D, 1);
output_args(3, :) = reshape(temp, [1, N]);
end

