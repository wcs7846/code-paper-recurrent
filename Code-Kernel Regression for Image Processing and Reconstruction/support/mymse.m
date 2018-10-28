function y = mymse(x)

y = sum(x(:).^2) / length(x(:));