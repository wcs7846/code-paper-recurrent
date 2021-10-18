function [ output_args ] = createEvidence( fe_array, lib )
%% CREATEEVIDENCE: the function of creating an evidence
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  fe_array    --> the array of the focal elements
         lib         --> a library of all possible options
 Output: output_args --> the evidence
 [Tips]: "type" only have 6 classes:
          (1) T: Triangle; (2) S: Square; (3) R: rectangle; (4) B: Bracket;
          (5) C: circle;   (6) D: disk;
         "curvature" only have 2 classes:
          (1) P: planes;   (2) C: curved;
%}
[~, N] = size(fe_array);

% [Tips]: the F represent the full set
% compute the BPA of F 
tmpFE.basicElement = lib;
tmpFE.bpa = 0;
for n = 1:N
    tmpFE.bpa = tmpFE.bpa + fe_array(n).bpa;
end
tmpFE.bpa = 1 - tmpFE.bpa;

% integrate all focal element to form the evidence
output_args = [fe_array, tmpFE];
end
