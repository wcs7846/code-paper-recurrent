function [ output_lsk ] = LSK( coordinate, row, col, gradient_x, gradient_y )
% LSK: local steering kernel algorithm
% reference: Robust infrared small target detection using local steering kernel reconstruction
% Copyright:2018-9-3 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
%{
 detail
 model: please refer the orginal article.
 Input:  coordinate --> the coordinate matrix of Image
         row        --> the rows of Image
         col        --> the cols of Image
         gradient_x --> the gradient matrix of patch(x axis)
         gradient_y --> the gradient matrix of patch(y axis)
 Output: output_lsk --> the LSK of patch
 Tips: the rows and cols of coordinate must be odd number
%}
%% Body: calculate the LSK of the patch
% step 0:set some constant parameters
lamda1 = 1; lamda2 = 10^(-7); alf = 3; h = 0.2;
[len,~] = size(coordinate)  ; s = sqrt(len);
padding = (s-1)/2;
K = zeros(row, col);
% step 1:calculate the matrix G: the first derivatives along the horizontal and vertical axes
Gx_vector = reshape(gradient_x, [], 1);
Gy_vector = reshape(gradient_y, [], 1);
G = [Gx_vector, Gy_vector];
% step 2:calculate the matrix C
% -- use SVD to get u1,u2 and s1,s2
[~,S,V] = svd(G, 'econ'); % Tips£º U and V is unitary matrix; S is diagonal matrix
s1 = S(1,1); s2 = S(2,2);
% Tips: (SVD)G = USV'; so the singular vector is V's row vector
u1 = V(1,:)'; u2 = V(2,:)';
% Tips: lamda1 = 1; lamda2 = 10^(-7); alf = 3({3,5,7} is available)
N = len;
t1 = ((s1*s2+lamda2)/N)^alf; % front part
t2 = ((s1+lamda1)/(s2+lamda1))*(u1*u1') + ((s2+lamda1)/(s1+lamda1))*(u2*u2'); % latter part
C = t1*t2;
% step 3:calculate the K
coeff = sqrt(det(C))/(h^2); % coeff = coefficient of K
% -- construct the delta_x vector: delta_x = x_j - x_i(x_i is central point)
% ncol: x coordinates; nrow: y coordinates;
% [ pVector, center_point, len ] = pointVector(ncol,ncol+s-1,nrow,nrow+s-1);
for n = 1:len
    tmp = double(coordinate(n,:))*C*double(coordinate(n,:))';
    tmp_K = exp(tmp/(-2*h*h));
    % fill the value of K into matrix K: pVector(n,:) is nth point's coordinates
    K(coordinate(n,1)+padding+1,coordinate(n,2)+padding+1) = coeff*tmp_K;
end
% step 4:output the result
output_lsk = K;
end

