function [ output_args ] = proximal_log( input_args, tao, sita, lamda )
%% PROXIMAL_LOG 
% reference: Unsupervised ridge detection using second order anisotropic Gaussian kernels
% Copyright:2019-2-27 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
% this is a auxiliary function to implement the proximal operator algorithm
%{
 detail
 Input:  input_args  --> the input image(2D/1D, and double type);
         tao         --> the normalization paratmer
         sita        --> the parameter in the log(|sigma|+sita)
         lamda       --> the scale factor of log function(default 1)
 Output: output_args --> the optimal(2D/1D, double)
%}
%% get the basic information
N_sig = length(input_args);
sig_X = zeros(size(input_args));
sig_Y = input_args;

%% calcualte some parameter;
new_sig_Y =  sig_Y - sita;
new_tao = lamda*tao - sig_Y*sita;

%% main loop
for n = 1:N_sig
    % judge
    t_judge = new_sig_Y(n)^2 - 4*new_tao(n);
    
    if t_judge >= 0
        % part 2
        t_sig = new_sig_Y(n) - 2*new_tao(n)/(new_sig_Y(n)+sqrt(t_judge));
        sig_X(n) = max(t_sig, 0);
    else
        % part 1
        sig_X(n) = 0;
    end
end
%% output 
output_args = sig_X;
end

