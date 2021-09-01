function [ cf, rmsf, ibw] = centerFreq( tfs, f )
%CENTERFREQ extract the center frequence from the time-frequence spectrum
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  tfs          --> the time-frequence spectrum (dim: 2D)
         f            --> the frequence vector (dim: 1D)
 Output: cf   --> the center frequence
         rmsf --> the root mean square frequence
         ibw  --> instantaneous bandwidth
 Ref:王雨青. 地震信号分数域属性分析及储层流体识别关键技术研究[D].电子科技大学,2018.
%}
p = abs(tfs).^2;
[~, M] = size(p); % M is number of time points
f_padd = repmat(f', 1, M);

numerator = sum(p.*f_padd, 1);
denumerator = sum(p, 1);
cf = numerator./denumerator; 

numerator = sum(p.*f_padd.*f_padd, 1);
rmsf = sqrt(numerator./denumerator);

f_diff = f_padd - cf(1);
numerator = sum(p.*f_diff.*f_diff, 1);
ibw = sqrt(numerator./denumerator);
end

