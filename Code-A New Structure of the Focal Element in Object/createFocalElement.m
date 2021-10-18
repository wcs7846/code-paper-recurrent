function [ output_args ] = createFocalElement( be_array, prop )
%% CREATEFORCEELEMENT: the function of creating a focus element
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  be_array    --> the array of basic elements
         prob        --> the BPA of this focal element
 Output: output_args --> the focal element
 [Tips]: "type" only have 6 classes:
          (1) T: Triangle; (2) S: Square; (3) R: rectangle; (4) B: Bracket;
          (5) C: circle;   (6) D: disk;
         "curvature" only have 2 classes:
          (1) P: planes;   (2) C: curved;
 [Example]:
 pf1 = createFocalElement(p_4tp, 0.5);
 pf1 = createFocalElement([p_4tp, p_1sp], 0.5);
%}
output_args.basicElement = be_array;
output_args.bpa = prop;
end

