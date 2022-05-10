function [mt,mtn,yprev,z] = tbmreg_val(x,y,w,xt,K,alpha,g,ymin,ymax);

% Propagation algorithm for TBM regression method
% [mt,mtn,yprev,z] = tbmreg_val(x,y,w,xt,K,alpha,g,ymin,ymax);

% Inputs:

% x: matrix of input vectors, size (n,d), n=nb pf examples, d=nb of features
% y: vector (n,1) of target outputs
% w: vector (n,1) or reliabilities of target output values (1: the value is fully reliable; 0: the value is actually
%    totally unknown)
% xt: matrix (nt,d) of test input vectors
% K: number of neighbors
% alpha: parameter used to compute the bba's (usually >= 0.9)
% g: gamma parameter (a scalar)
% ymin, ymax: lower and upper bounds for the domain of the output variable

% Outputs:

% mt: matrix (nt,c+1) of unnormalized output bba's for teh test examples
% mtn: matrix (nt,c+1) of normalized output bba's for teh test examples
% yprev: vector (nt,1) of point predictions for the test examples
% z: vector of the c distinct values of y in the training sets (the focal elements of the output bba's)

% USES:
% err_reg.m : error function, returns the error and the gradient
% harris.m : optimization procedure

% Reference :
% S. Petit-Renaud and T. Denoeux. Nonparametric regression analysis of uncertain and imprecise data using 
% belief Functions. International Journal of Approximate Reasoning, Vol. 35, No. 1, 1-28, 2004.

[n,d]=size(x); 
nt=size(xt,1);

[z,ii,jj]=unique(y);
c=length(z);
m=zeros(n,c+1);
for i=1:n,
    m(i,jj(i))=w(i);
    m(i,c+1)=1-w(i);
end;
z(c+1)=(ymax+ymin)/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the K-nearest neighbors in the training set

D=xt*x';
D=repmat(sum(xt.^2,2),1,n)+repmat(sum(x.^2,2)',nt,1)-2*D;
R=diag(Inf*ones(n,1));

ist=zeros(nt,K);
for i=1:K
    [tmp,I]=min(D,[],2);
    ist(:,i)=I;
    D=max(D,R(I,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the bba's

mk = [zeros(nt,c) ones(nt,1)];
mt = mk;
phi = zeros(nt,K);

for k = 1:K
 	I = ist(:,k);
    phi = alpha*exp(-g^2*sum((x(I,:)-xt).^2,2));
 	mk = [(m(I,1:c).*repmat(phi,1,c)) (1-(1-m(I,c+1)) .* phi)];
	mt = [(mt(:,1:c).*(mk(:,1:c)+repmat(mk(:,c+1),1,c))+...
            mk(:,1:c).*repmat(mt(:,c+1),1,c)) (mt(:,c+1).*mk(:,c+1))];
end

% normalization

Kn = sum(mt,2);
mtn = mt./(Kn*ones(1,c+1));
yprev=mtn*z;




