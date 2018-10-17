function [ output ] = LHM( Input_Img, s )
% LHM: local homogeneity measure algorithm
% reference: An infrared small target detection method based on multiscale local homogeneity measure
% operation
% Copyright:2018-9-3 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
%{ 
 detail
 model: please refer the orginal article.
 Input:  Input_Img --> the matrix of Image(must be 2D)
         s         --> the size of central patch(must be odd number)
 Output: output    --> the feature map of local homogeneity measure
 
 Tips: the input matrix must be 2D matrix
%}

% step 1: get the patch image
[row, col] = size(Input_Img);
boundary = (s-1)/2+s;
matrix_sita = zeros(row, col);
matrix_w    = zeros(row, col);
for nrow = 1:row-3*s
    for ncol = 1:col-3*s   
        % patch image
        patch = double(Input_Img(nrow:nrow+3*s-1,ncol:ncol+3*s-1));
        % calculate the standard deviation of central patch
        central_patch = patch(1+s:2*s,1+s:2*s);
        sita_row = boundary + nrow;
        sita_col = boundary + ncol;
        sita = std(central_patch(:));
        matrix_sita(sita_row, sita_col) = sita;
        % calculate the Inter-patch heterogeneity measure
        meanVector = zeros(1,9);
        % calculate the D vectore and use the block matrix
        for t = 1:9
            % Tips: The reason to set the tt and tp temporal parameter is at
            % the bottom
            tt = floor((t-1)/3);
            tp = mod(t-1,3);
            temp_block = patch(1+tt*s:(tt+1)*s, 1+tp*s:(tp+1)*s);
            meanVector(t) = mean(temp_block(:));
        end
        D_temp = [meanVector(1:4), meanVector(6:9)];
        % calculate the mean of the central patch (the calculate order can be changed)
        mT = meanVector(5);
        %
        D = -(D_temp - mT);
        % calculate the value C
        di = zeros(1,4);
        for t = 1:4
            di(t) = D(t)*D(t+4);
            if di(t) < 0
                di(t) = 0;
            end
        end
        C = min(di);
        % calculate the value w
        matrix_w(sita_row, sita_col) = C*(sita+0.01)^-1;
    end
end

% step 2: output the result
output = matrix_w(boundary:boundary+row-3*s, boundary:boundary+col-3*s);
end
%% the explaination for why set temporal parameter: tt and tp
% an analysis of 3 orders
% t  = 1,2,3,4,5,6,7,8,9
% tt = 0,0,0,1,1,1,2,2,2
% tp = 0,1,2,0,1,2,0,1,2
% So: tt = floor((t-1)/3); tp = mod((t-1),3);
