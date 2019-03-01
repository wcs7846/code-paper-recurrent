function [e, gMag] = nonmaxSuppress(varargin)
%NONMAXSUPPRESS does non-maximum suppression on 2-D and 3-D arrays
%   E = NONMAXSUPPRESS(G, NAME, VALUE, ...) finds directional local maxima
%   in a gradient array, as used in the Canny edge detector, but made
%   separately accessible here for greater flexibility. The result is a
%   logical array with the value true where the gradient magnitude is a
%   local maximum along the gradient direction.
%
%   G is a cell array containing the components of the gradient, of length
%   2 or 3. All the elements of G must have the same size; they must be 2-D
%   if the length of G is 2 or 3-D if the length of G is 3.
%
%   [E, GMAG] = NONMAXSUPPRESS(G, NAME, VALUE, ...) also returns an array
%   containing the magnitude of the gradient at each point.
%
%   Name-value pairs allow control over the process. Parameter names are:
%
%   'Method' - Specifies the algorithm used for interpolation amongst
%   neighbouring elements of the gradient magnitude array. For a point to
%   be marked as an edge, the gradient magnitude must be greater than the
%   interpolated magnitudes at two nearby points, one along the gradient
%   direction and the other in the opposite direction. Methods can be:
%
%       'boundary' (default) - Linear interpolation is carried out on the
%       boundary of the 2x2 square (in 2-D) or the 2x2x2 cube (in 3-D)
%       centred on the point being tested.
%
%       'nearest', 'linear', 'cubic', 'spline' - Interpolation is carried
%       out in 2-D or in 3-D. The method is passed to GRIDDEDINTERPOLANT.
%
%   'Radius' (default 1) - The distance along the gradient direction at
%   which the interpolation for each point is done. A value close to 1 is
%   recommended. This parameter is ignored if the 'boundary' method is
%   specified.
% 
%   'SubPixel' (default false) - If true, sub-pixel interpolation is
%   carried out on the edge positions by fitting a parabola to the three
%   gradient magnitues round the peak. The result E will be a structure
%   with two fields:
%
%       E.edge - a logical array with true at the edge positions
% 
%       E.subpix - In the 2-D case, if E.edge(I,J) is true,
%       E.subpix{1}(I,J) is the first (row) coordinate of the estimated
%       sub-pixel position of the edge point. E.subpix{2} gives the second
%       (column) coordinate. The 3-D case is similar.
% 
% See also: canny, griddedInterpolant
% Argument sorting out
inp = inputParser;
inp.addRequired('g', ...
    @(g) validateattributes(g, {'cell'}, {'vector'}));
inp.addParameter('Method', 'boundary', ...
    @(s) ismember(s, {'boundary' 'linear' 'nearest' 'cubic' 'spline'}));
inp.addParameter('Radius', 1, ...
    @(s) validateattributes(s, {'numeric'}, ...
    {'real' 'positive', 'nonnan', 'finite', 'scalar'}));
inp.addParameter('SubPixel', false, ...
    @(s) validateattributes(s, {'logical'}, {'scalar'}));
inp.parse(varargin{:});
g = inp.Results.g;
method = inp.Results.Method;
rad = inp.Results.Radius;
% Get and check no. dimensions. Currently only to 2 or 3-D, but could
% generalise to N if needed.
ndim = length(g);
if ndim < 2 || ndim > 3
    error('DavidYoung:nonmaxSuppression:dimsNot2or3', ...
        'Length of gradients array must be 2 or 3, was %d', ndim);
end
% make the gradients into an array
gcomps = cell2cols(g);
% compute gradient magnitude array
gMagCol = sqrt(sum(gcomps.^2, 2));
gMag = reshape(gMagCol, size(g{1}));
% convert gradients to offsets
if strcmp(method, 'boundary')
    % Make largest offset 1 - on boundary of square or cube.
    % Might expect that specialised code would be faster than general
    % interpolation in this case. Not so.
    [~, maxgradcomp] = max(abs(gcomps), [], 2);
    gnorm = gcomps(...
        sub2ind(size(gcomps), (1:size(gcomps,1)).', maxgradcomp));
    method = 'linear';   % linear interp on square or cube boundary
else
    % put on sphere of given radius
    gnorm = gMagCol / rad;
end
gcomps = bsxfun(@rdivide, gcomps, gnorm);
% coords of data points in same format
coords = arrayfun(@(s) 1:s, size(gMag), 'UniformOutput', false);
[coords{:}] = ndgrid(coords{:});
coords = cell2cols(coords);
% get interpolant
interpolant = griddedInterpolant(gMag, method);
% interpolate in gradient direction and its opposite
ginta = interpolant(coords + gcomps);
gintb = interpolant(coords - gcomps);
% suppress nonmax pixels
ok = gMag(:) >= ginta & gMag(:) >= gintb;
edges = reshape(ok, size(gMag));
if inp.Results.SubPixel
    e.edge = edges;
    e.subpix = subpix(gMag, ginta, gintb, ok, coords, gcomps);
else
    e = edges;
end
end
function a = cell2cols(c)
% put the arrays in cell array c into the columns of a
a = cell2mat(cellfun(@(x) x(:), c, 'UniformOutput', false));
end
function pos = subpix(gMag, ginta, gintb, ok, coords, gcomps)
% points on parabola
q = gMag(ok);
r = ginta(ok) - q;    % point in direction of gradient
p = gintb(ok) - q;    % opposite direction to gradient
offset = (p-r)./(2*(p+r)); % scalar offset of peak 
offset = bsxfun(@times, offset, gcomps(ok, :)); % vector offset of peak
offcoords = coords(ok, :) + offset;
res = zeros(size(gMag));
for i = 1:size(coords,2)
    res(ok) = offcoords(:, i);
    pos{i} = res;
end
end