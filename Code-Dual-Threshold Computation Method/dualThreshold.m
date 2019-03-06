function [ upper_thres, lower_thres ] = dualThreshold( input_args, scale )
%% DUALTHRESHOLD 
% reference: Bowstring-Based Dual-Threshold Computation Method for Adaptive Canny Edge Detector
% Copyright:2019-2-28 MarkLHF, UESTC.(e-mail:2751867750@qq.com)
%{
 detail
 Input:  input_args  --> the gradient magnitude image(2D image)
         scale       --> the scale of data type 
                         (Example: double = [0,1] and scale = 1)
 Output: upper_thres --> the upper threshold of image(for binary) 
         lower_thres --> the lower threshold of image
%}
%% <1> calculate the gradient magnitude(GMI)
gmi = input_args; % due to this gradient I had calculated before

%% <2> normalize the GMI
% [Tips]: GMI(x,y) >= 0
norm_coeff_gmi = max(gmi(:));% the coefficient of normalization (GMI)
gmi_norm = (gmi./norm_coeff_gmi).*scale;

%% <3> build and normalize the GMH 
gmh = imhist(gmi_norm);

norm_coeff_gmh = max(gmh(:)); % the coefficient of normalization (GMH)
gmh = (gmh./norm_coeff_gmh).*scale;

%% <4.1> determine the slope of the long-/short-bowstring.
% search the single peak point 
peakValue = max(gmh);
peakPoint = find(gmh == peakValue);

% determine almost zero point(AZP)
gamma = 0.005 * 0.01;
azpValue = gamma * peakValue;
azp = find(gmh <= azpValue);
azp = azp(length(azp)); % set the rightmost point as AZP 
% determine kH parameter
kH = (gmh(peakPoint) - gmh(azp))/(peakPoint - azp);
%% <4.2> obtain the temporary upper-threshold 
i_range = peakPoint:length(gmh); 
lH = gmh(i_range)' - kH.*i_range;

iH = find(lH == max(lH));
iH = peakPoint + iH - 1;
%% <4.3> obtain the temporary lower-threshold
% determine kL parameter
kL = (gmh(peakPoint) - gmh(iH))/(peakPoint - iH);

i_range = peakPoint:2:iH; 
lL = gmh(i_range)' - kL.*i_range;

iL= find(lL == max(lL));
iL = peakPoint + iL - 1;
%% <5> refine the dual-threshold
% I have some problem in this part,  
K = length(gmh);
% output the result
upper_thres = iH/K*scale;
lower_thres = iL/K*scale;
end

