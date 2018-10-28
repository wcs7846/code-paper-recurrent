function y = round0_255(x)

y = x;
below0 = y < 0;
above255 = y > 255;

y = y .* abs(below0 - 1);
y = y .* abs(above255 - 1) + above255 * 255;