function [ y_pred ] = lssvr_predict( x_input, x_train, alfa, bias, g )
%LSSVR_PREDICT -- instruction
%   details:
%[Input]
%   x_input: the input range x; like x in the f(x)
%   x_train: the train x
%   alfa: Lagrange multipliers, N*1 vector, N is the number of sample
%   bias: a bias of the linear model, a number
%   g: Radial Basis Functon learning parameter, is equal to 1/2sigma^2(gaussian function)
%[Output]
%   y_pred: the prediected value
%[Reference]
%   https://www.mathworks.com/matlabcentral/fileexchange/73706-support-vector-regression-ls-svr-for-non-linear-functions
%   [1]	L. T. Qin, S. S. Liu, H. L. Liu, and Y. H. Zhang, "Support vector regression and least squares support vector regression for hormetic dose-response curves fitting," (in English), Chemosphere, vol. 78, no. 3, pp. 327-334, Jan 2010, doi: 10.1016/j.chemosphere.2009.10.029.
%   Ludovico Cuoghi (2022). Support Vector Regression ( LS-SVR) for Non Linear functions (https://www.mathworks.com/matlabcentral/fileexchange/73706-support-vector-regression-ls-svr-for-non-linear-functions), MATLAB Central File Exchange. Retrieved April 22, 2022.
if nargin < 4
    g = 0.01;
end
[~,x_length] = size(x_input);
xs_length=length(x_train);

A = x_input'*x_train;
M  = diag(diag(x_input'*x_input))*ones(x_length,xs_length);
M2 = ones(x_length,xs_length)*diag(diag(x_train'*x_train));
A=M+M2-2*A;
A = exp(-g*A);

y_pred = A*alfa+bias;
end

