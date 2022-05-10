function [ alfa, bias ] = lssvr_train( x_train, y_train, C, g )
%LSSVR_TRAIN -- instruction
%   details:
%[Input]
%   x_train: the input range x; like x in the f(x); each x should be a col vector
%   y_train: the input data
%   C: Parameter defined to avoid overfitting
%   g: Radial Basis Functon learning parameter, is equal to 1/2sigma^2(gaussian function)
%[Output]
%   alfa: Lagrange multipliers, N*1 vector, N is the number of sample
%   bias: a bias of the linear model, a number
%[Reference]
%   https://www.mathworks.com/matlabcentral/fileexchange/73706-support-vector-regression-ls-svr-for-non-linear-functions
%   [1]	L. T. Qin, S. S. Liu, H. L. Liu, and Y. H. Zhang, "Support vector regression and least squares support vector regression for hormetic dose-response curves fitting," (in English), Chemosphere, vol. 78, no. 3, pp. 327-334, Jan 2010, doi: 10.1016/j.chemosphere.2009.10.029.
%   Ludovico Cuoghi (2022). Support Vector Regression ( LS-SVR) for Non Linear functions (https://www.mathworks.com/matlabcentral/fileexchange/73706-support-vector-regression-ls-svr-for-non-linear-functions), MATLAB Central File Exchange. Retrieved April 22, 2022.
%% parameter check
if nargin < 2
    C = 100;
elseif nargin < 3
    g = 0.01;
end
%% main body
xs_length=length(x_train);
A = x_train'*x_train;
M=diag(diag(A))*ones(xs_length,xs_length); % -> N = repmat(diag(A2),[1,N]) 
A=M+M'-2*A;
% R=diag(Inf*ones(xs_lenght,1)); % remove the self
% A2=max(A2,R); % A2 means the distrances of these x (||x-xi||^2)
A = exp(-g*A);
A = A + 1/C*diag(ones(xs_length,1));

O=[0, ones(1,xs_length);
        ones(xs_length,1), A];

b=zeros(xs_length+1,1);
c=zeros(xs_length+1,1);

c(2:xs_length+1) = y_train';

b=inv(O)*c; %it contains LaGrange multipliers and bias

% For convenience, we separate the LaGrange Multipliers and the Bias term

%Bias term
bias=b(1,1);

%Define LaGrange Multipliers alfa
alfa=b(2:xs_length+1); 

end

