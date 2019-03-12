function [ output_X, output_E, list_e_norm ] = arllre( input_args, itera, show_debug, lamda)
%% ARLLRE
% reference: Unsupervised ridge detection using second order anisotropic Gaussian kernels
% Copyright:2019-2-27 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
%{
 detail
 Input:  input_args --> the input image(2D, and double type);
         itera      --> the number of iteration
         show_debug --> the switch to show the debug information
 Output: output_X   --> the reconstruction(2D, double)
         output_E   --> the residual
         list_e_norm --> the list of the norm of E matrix
%}
%% get the basic information and allocate the store space
[row, col] = size(input_args);

sita = 10^-4; 
% sigma = 1/(5*std(input_args(:)));
% tao = 1.2*sqrt(2)*sigma^2;
mu = 1 / (5*std(input_args(:)));

X = zeros(row, col);
W = zeros(row, col);
E = zeros(row, col);

list_e_norm = zeros(itera, 1);
%% main loop
for k = 1:itera
    tic;
    %% update the X_k+1 matrix
    % get matrix of Z_X
    Z_X = input_args - E + W./mu;
    % get the singular value matrix of Z_k+1
    [U, t_sig_Y, V] = svd(Z_X);
    sig_Y = diag(t_sig_Y);

%     sigma = 1/(5*std(Z_X(:)));
    tao = 1/mu;
    % calculate the every singular by proximal operator algorithm
    sig_X = proximal_log(sig_Y, tao, sita, 1);
    
    % reconstruct
    S_x = diag(sig_X);   
    if row >= col
        d = row - col;
        S_x = [S_x; zeros(d, col)];
    else
        d = col - row;
        S_x = [S_x, zeors(row, d)];
    end
    newX = U * S_x * V';
%     newX = singularValueShrinkage(Z_X, tao);
    %% update the E matrix 
    % get matrix of Z_E
    Z_E = input_args - newX + W./mu;
    % calculate the every singular by proximal operator algorithm
    %     Z_E = abs(Z_E) + sita;
    %     E = softThreshold(Z_E .* (Z_E>0), sita);
%         lamda = 100000;  
    
    t_ze = reshape(Z_E, [row*col, 1]);
    t_E = proximal_log(t_ze .* (t_ze>0), tao/lamda, sita, 1); % tips: sig_E
%     newE = softThreshold(Z_E .* (Z_E>0), lamda/mu);
    % reconstruct
    newE = reshape(t_E, [row, col]);
    %% update the W matrix and mu
    newW = W + mu*(input_args - newX - newE);
    new_mu = 1.5*mu;
    
    elapsedTime = toc;
    %% record and stop criterion
    stopCriterion = norm(input_args - newX - newE, 'fro') / norm(input_args, 'fro');
    list_e_norm(k) = stopCriterion;
    if show_debug
        disp(sprintf("[#Iter %2d] time is %f, ||E||_f: %f", k, elapsedTime, list_e_norm(k)));
        pause(0.01);
    end
    if stopCriterion < 1e-5
        break;
    end
    %% Assignment to continue
    X = newX;    E = newE;
    W = newW;    mu = new_mu;
end
%% output
output_X = X;
output_E = E;
end

