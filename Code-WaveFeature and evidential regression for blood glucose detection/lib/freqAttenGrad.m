function [ fag ] = freqAttenGrad( tfs, fvector )
%% FreqAttenGrad: Frequency attenuation gradient calculation function
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  tfs     --> the time-frequence spectrum (dim: 2D)
         fvector --> the frequecy vector
 Output: fag     --> the frequency attenuation gradient
 Ref:王雨青. 地震信号分数域属性分析及储层流体识别关键技术研究[D].电子科技大学,2018.
%}
% Tips: the row represent every frequency point
%       the col represent every time point
%% body
energy = sum(abs(tfs), 1);
[Nfreq, len] = size(tfs);
fag = zeros(1, len);
% calculate the 0.65 and 0.85 frequency point
for n = 1 : len
    tol_energy = energy(n); % the total energy at this time point
    tmp_energy = abs(tfs(:, n));
    cumu_energy = zeros(1, Nfreq); tt = 0;
    for m = 1:Nfreq
        tt = tt + tmp_energy(m);
        cumu_energy(m) = tt; % Cumulative energy distribution
    end
    % 0.65 frequency point --> M
    diff = abs(cumu_energy - 0.65*tol_energy);
    [~, minIndex1] = min(diff);
    % 0.85 frequency point    
    diff = abs(cumu_energy - 0.85*tol_energy);
    [~, minIndex2] = min(diff);  
    
    range = min(minIndex1, minIndex2):max(minIndex1, minIndex2);
    xdata1 = fvector(range);
    ydata1 = abs(tfs(range, n ));
    
    p1 = polyfit(xdata1,ydata1',1);
    fag(n) = p1(1);
end
fag(fag>0) = 0;
fag = fag/max(abs(fag));
end

