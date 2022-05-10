function [g,err,z] = tbmreg_fit(x,y,w,K,alpha,g0,ymin,ymax);

% Learning algorithm for TBM regression method
% [g,err,z] = tbmreg_fit(x,y,w,K,alpha,g0,ymin,ymax);

% Inputs:

% x: matrix of input vectors, size (n,d), n=nb pf examples, d=nb of features
% y: vector (n,1) of target outputs
% w: vector (n,1) or reliabilities of target output values (1: the value is fully reliable; 0: the value is actually
%    totally unknown)
% K: number of neighbors
% alpha: parameter used to compute the bba's (usually >= 0.9)
% g0: initial guess for gamma parameter (a scalar)
% ymin, ymax: lower and upper bounds for the domain of the output variable

% Outputs:

% g: optimized scalar gamma parameter
% err: training error
% z: vector of the c distinct values of y in the training sets (the focal elements of the output bba's)

% USES:
% err_reg.m : error function, returns the error and the gradient
% MATLAB optimization toolbox (function fminbnd)

% Reference :
% S. Petit-Renaud and T. Denoeux. Nonparametric regression analysis of uncertain and imprecise data using 
% belief Functions. International Journal of Approximate Reasoning, Vol. 35, No. 1, 1-28, 2004.

[n,d]=size(x);
[z,ii,jj]=unique(y);
c=length(z);
m=zeros(n,c+1);
for i=1:n
    m(i,jj(i))=w(i);
    m(i,c+1)=1-w(i);
end
z(c+1)=(ymax+ymin)/2;

% Search for K nearest neighbors in training set

D=x*x';
N=diag(diag(D))*ones(n,n); % -> N = repmat(diag(D),[1,N]) 
D=N+N'-2*D;
R=diag(Inf*ones(n,1)); % remove the self
D=max(D,R); % D means the distrances of these x (||x-xi||^2)

is=zeros(n,K); % the minimum in the first col.
for i=1:K
    [tmp,I]=min(D,[],2); % the minimum of each row.
    is(:,i)=I;
    D=max(D,R(I,:));
end
% Optimization


options=optimset('fminbnd');
g=fminbnd('error_reg',eps,20,options,x,y,m,z,is,K,alpha);
err=error_reg(g,x,y,m,z,is,K,alpha);







