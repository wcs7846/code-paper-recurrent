function x = Ninv(p)
% inverse normal
x = -sqrt(2).*erfcinv(2*p);
end