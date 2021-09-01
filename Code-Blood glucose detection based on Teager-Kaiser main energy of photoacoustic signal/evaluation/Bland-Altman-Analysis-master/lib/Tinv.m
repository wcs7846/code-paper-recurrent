function x = Tinv(p,n)
% inverse t
assert(p>=.75 & p<=1,'p must satisfy 0.75<=p<=1.')
if n<=1
    error 'Not implemented for n<=1.'
elseif n<1000
    bii = betaincinv(2*abs(p-1/2), n/2, 1/2, 'upper');
    x = sqrt(n*(1-bii)/bii);
else % n>=1000
    % Reference:
    % Abramowitz & Stegun, 10th ed. 1972, formula 26.7.5
    xp = Ninv(p);
    g1 = @(x) (x^3+x)/4;
    g2 = @(x) (5*x^5+16*x^3+3*x)/96;
    g3 = @(x) (3*x^7+19*x^5+17*x^3-15*x)/384;
    g4 = @(x) (79*x^9+776*x^7+1482*x^5-1920*x^3-945*x)/92160;
    x = xp + g1(xp)/n + g2(xp)/n^2 + g3(xp)/n^3 + g4(xp)/n^4;
end
end