function [st,t,f] = st(x,p)
minfreq = 0;
maxfreq = fix(length(x)/2);  %向0取整
t = (0:length(x)-1);  %时间间隔取1
fre_number =ceil((maxfreq - minfreq+1));  %向上取整，频率间隔取1
f = (minfreq + [0:fre_number-1])/(length(x));% 归一化的频率，显示时乘以采样频率即可与真实频率一一对应
st = strans(x,minfreq,maxfreq,p); 
return

function st = strans(x,minfreq,maxfreq,p)
n=length(x);
vector_fft=fft(x);
vector_fft=[vector_fft,vector_fft];%在平移时也包含完整的一个信号
st=zeros(ceil((maxfreq - minfreq+1)),n);
st(1,:) = mean(x)*(1&(1:n));  % Tips: this is  the average of signal when frequency is zero
for i=1:(maxfreq-minfreq)
   st(i+1,:)=ifft(vector_fft(minfreq+i+1:minfreq+i+n).*g_window(n,minfreq+i,p));
end   

% function gauss=g_window(length,freq)
% vector(1,:)=[0:length-1];
% vector(2,:)=[-length:-1];
% vector=vector.^2;    
% vector=vector*(-2*pi^2/freq^2);
% gauss=sum(exp(vector));
function gauss=g_window(length,freq,p)
vector(1,:)=[0:length-1]; % the left  side
vector(2,:)=[-length:-1]; % the right side 
vector=vector.^2;    
vector=vector*(-2*pi^2*p^2/freq^2);
% vector=vector*(-2*pi^2/freq^(2*p));
gauss=sum(exp(vector));