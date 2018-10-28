function C = steering(zx, zy, I, wsize, lambda, alpha)
% Compute steering matrices
%
% [USAGE]
% C = steering(zx, zy, I, wsize, lambda, alpha)
%
% [RETURNS]
% C      : steering matrices
%
% [PARAMETERS]
% zx, zy : image gradients along x and y directions
% I      : sampling positions
% wsize  : size of an analysis window
% lambda : regularization parameter
% alpha  : structure sensitive parameter
%
% [HISTORY]
% Nov  19, 2005 : created by hiro
% June 16, 2005 : modified by hiro
% Apr  14, 2008 : economical SVD, by hiro

[N, M] = size(zx);
C = zeros(2, 2, N, M);

if mod(wsize, 2) == 0
    wsize = wsize + 1;
end
win = floor(wsize / 2);

K = fspecial('disk', win);
K = K ./ K(win+1, win+1);
% K = ones(wsize,wsize);

% mirroring
zx = EdgeMirror(zx, [win, win]);
zy = EdgeMirror(zy, [win, win]);

for i = 1 : N
    for j = 1 : M
        
        if I(i,j) == 0
            continue;
        end
        
        gx = zx(i:i+wsize-1, j:j+wsize-1) .* K;
        gy = zy(i:i+wsize-1, j:j+wsize-1) .* K;

        G = [gx(:), gy(:)];
        len = sum(K(:));
        [u s v] = svd(G, 0);
        S(1) = (s(1,1) + lambda) / (s(2,2) + lambda);
        S(2) = (s(2,2) + lambda) / (s(1,1) + lambda);
        tmp = (S(1) * v(:,1) * v(:,1)' + S(2) * v(:,2) * v(:,2)')  * ((s(1,1) * s(2,2) + 0.0000001) / len)^alpha;
        C(1:2, 1:2, i, j) = tmp;
    end
end