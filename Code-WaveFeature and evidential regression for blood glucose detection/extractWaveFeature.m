function [ output_args ] = extractWaveFeature( x )
%EXTRACTWAVEFEATURE 此处显示有关此函数的摘要
%   此处显示详细说明
%INPUT
% x:
%OUTPUT
% output_args: [wn]
[NumSample,len]=size(x);

meanValue = mean(x,2); stdValue = std(x,0,2);
[maxValue, loc] = max(abs(x-repmat(meanValue,[1,len])),[],2);

Ft = normpdf(1:len, round((loc+loc))/2, 1);Ft = Ft./repmat(max(Ft, [], 2), [1, len]); %log10(len)-1
Ft = Ft.*repmat(maxValue, [1, len])+repmat(meanValue, [1,len]);
% basic equation parameter
wn = zeros(NumSample, 1);
options=optimset('fminbnd');func = @error_reg;
for n = 1:NumSample
    data = x(n,:); xrange = (0:len-1);
    wn_opt=fminbnd(func, 0,1, options, xrange, data, Ft(n,:));
    wn(n) = wn_opt;   
    % show
%     rho = 0.9;
%     figure;plot(xrange, data, 'r-');hold on;
%     [t2,y] = ode45(@(t,y)vdp1(t,y,wn_opt,rho,Ft(n,:)),[min(xrange) max(xrange)],[mean(data);0]);
%     yi= interp1(t2,y,1:len,'spline');
%     plot(xrange, yi(:,1), 'b-');hold on;
end
% Maximum overshoot
Mp = (maxValue - meanValue);
% Adjustment time
ts = 3./wn;

output_args = [wn, Mp, ts];
end

function err = error_reg(wn, t, testdata, Ft)
rho = 0.9;len = length(testdata); meanValue = mean(testdata);

[t2,y] = ode45(@(t,y)vdp1(t,y,wn,rho,Ft),[0 len-1],[meanValue;0]);
yi= interp1(t2,y,t,'spline'); % 'method'是最邻近插值， 'linear'线性插值； 'spline'三次样条插值； 'cubic'立方插值

[yupper,~] = envelope(testdata,50); [yupperi,~] = envelope(yi(:,1)',50);
% figure;
% subplot(211);plot(t, yi(:,1), 'r-');hold on;legend('信号波形');title(strcat('wn=',num2str(wn)));
% subplot(212);plot(t, yi(:,2), 'b-');hold on;legend('信号1阶导数');
err = norm(yupper-yupperi);
end

function dydt = vdp1(t,y,rho, wn, Ft)
dydt = [y(2); 
		-2*rho*wn*y(2)-wn^2*y(1)+Ft(round(t+1))];
end