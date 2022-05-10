function [err,G]=error_reg(g,x,y,m,z,is,K,alpha)
% Error function used by tbmreg_fit
% [err,G]=error_fun(g,x,y,m,z,is,K,alpha)
%
% Inputs:
% g: gamma parameter, to be optimized
% x: matrix of input vectors, size (n,d), n=nb pf examples, d=nb of features
% m: matrix of bbas, size (n,c+1); focal element c+1 is output frame
% is: indices of the K nearest neighbors
%     is(i,k) = index of the K-th NN of example i
% K: nb of neighbors
% alpha: parameter used to compute the bba's (usually >= 0.9)
%
% Outputs
% err= value of objective function
% G: derivative of error function with respect to g (empty in the current version; the derivative is evaluated
%    numerically by the optimization procedure harris.m

[n,d]=size(x);
c=size(m,2)-1;


% Error calculation

mk = [zeros(n,c) ones(n,1)];
mm = mk;
phi = zeros(n,K);

for k = 1:K
 	I = is(:,k);
    phi(:,k) = alpha*exp(-g^2*sum((x(I,:)-x).^2,2));
 	mk = [(m(I,1:c).*repmat(phi(:,k),1,c)) (1- (1-m(I,c+1)) .* phi(:,k))];
	mm = [(mm(:,1:c).*(mk(:,1:c)+repmat(mk(:,c+1),1,c))+...
	mk(:,1:c).*repmat(mm(:,c+1),1,c)) (mm(:,c+1).*mk(:,c+1))];
end
Kn = sum(mm,2);
mmn = mm./(Kn*ones(1,c+1));

yprev=mmn*z;
err=mean((y-yprev).^2);
G=[];
