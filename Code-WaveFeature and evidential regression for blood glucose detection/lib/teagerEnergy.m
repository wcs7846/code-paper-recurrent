function [ tk ] = teagerEnergy( tfs )
%% TEAGERENERGY: TK energy calculation function
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  tfs  --> the time-frequence spectrum (dim: 2D)
 Output: tk   --> the teager energy 
 Ref:王雨青. 地震信号分数域属性分析及储层流体识别关键技术研究[D].电子科技大学,2018.
%}
% Tips: the row represent every frequency point
%       the col represent every time point
%% body
rp_tfs = real(tfs);  % the real part
ip_tfs = imag(tfs); % the image part 

[row, col] = size(tfs);
tk = zeros(row, col);
for ncol = 2:col-1
    tk_r = rp_tfs(:, ncol).^2 - rp_tfs(:, ncol-1).*rp_tfs(:, ncol+1);
    tk_i = ip_tfs(:, ncol).^2 - ip_tfs(:, ncol-1).*ip_tfs(:, ncol+1);

    tk(:, ncol) = tk_r + tk_i;
end
tk = max(tk,[],1);
end

