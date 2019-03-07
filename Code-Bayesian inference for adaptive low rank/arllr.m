function [ X, E ] = arllr( input_args, lamda, ro )
%% ARLLR 
% reference: Unsupervised ridge detection using second order anisotropic Gaussian kernels
% Copyright:2019-2-27 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
%{
 detail
 Input:  input_args --> the input image(2D, and double type);
         lamda      --> the initial value
         ro         --> the step(a constant)
 Output: X --> the reconstruction(2D, double)
         E --> the residual
%}
%% get the basic information and allocate the store space
sita = 10^-4; 
sigma = std(input_args(:));
tao = 1.2*sqrt(2)*sigma^2;
total_itera = 1;

[row, col] = size(input_args);

X = zeros(row, col);
E = zeros(row, col);
% get the singular value matrix of Y and X

[U, t_sig_Y, V] = svd(input_args); 
sig_Y = diag(t_sig_Y);

sig_X = zeros(size(sig_Y));
%% main loop
for k = 1:total_itera
    %
    new_sig_Y =  sig_Y - sita;
    new_tao = tao - sig_Y*sita;
    
    N_sig = length(sig_Y);
    % calculate the every singular
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
    
    % reconstruct
    X = diag(sig_X);
    
    if row >= col
        d = row - col;
        X = [X; zeros(d, col)];
    else
        d = col - row;
        X = [X, zeors(row, d)];
    end
    X = U*X*V;
    
    E = input_args - X;
end


end

