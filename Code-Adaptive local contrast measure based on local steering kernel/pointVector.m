function [ output_args, center_point, length ] = pointVector( x_start, x_end, y_start, y_end )
% Summary: this function construct a vector of point with the first col is
% x coordinates and the second col is y coordinates
% Copyright:2018-9-6 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
%{ 
 detail
 Input:  x_start --> the start point of x
         x_end   --> the end point of x
         y_start --> the start point of y
         y_end   --> the end point of y
 Output: output_args  --> the vector of point
         center_point --> the vector of central point
         length       --> the length of output_args
 Tips: the input of number must be non-negetive integer
%}

% step 1: calculate the total number of the point in the vector
num_x = abs(x_end - x_start) + 1; % x is 
num_y = abs(y_end - y_start) + 1;
num_total = num_x * num_y;
length = num_total;
output_args = zeros(num_total, 2);
% step 2: fill the matrix -> output_args
% -- output_args: like [x,y]; the first col is x coordinates and the second
% is y
vector_y = y_start : y_end;
vector_x = x_start : x_end;
vector_t = ones(num_y, 1);
for n = 1:num_x
    output_args(1+(n-1)*num_y:n*num_y, 1) = vector_t * vector_x(n);
    output_args(1+(n-1)*num_y:n*num_y, 2) = vector_y;
end
center_point = [round((x_start+x_end)/2), round((y_start+y_end)/2)];
end

