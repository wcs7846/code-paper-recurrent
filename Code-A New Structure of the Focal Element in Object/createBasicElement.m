function [ output_args ] = createBasicElement( num_planes, type, curvature )
%% CREATEELEMENT: the function of creating a basic element
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  num_planes  --> the number of this type planes (a integer)
         type        --> the type of planes 
         curvature   --> the curvature of planes
         prob        --> the BPA of this element
 Output: output_args --> the basic element
 [Tips]: "type" only have 6 classes:
          (1) T: Triangle; (2) S: Square; (3) R: rectangle; (4) B: Bracket;
          (5) C: circle;   (6) D: disk;
         "curvature" only have 2 classes:
          (1) P: planes;   (2) C: curved;
 [Tips]: this element is the one part of focal elements
%}
output_args.num_planes = num_planes;
output_args.type = upper(type);
output_args.curvature = upper(curvature);
end

