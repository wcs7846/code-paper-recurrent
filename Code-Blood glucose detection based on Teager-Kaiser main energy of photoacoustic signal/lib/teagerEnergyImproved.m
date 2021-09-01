function [ tk ] = teagerEnergyImproved( tfs )
%% TEAGERENERGY: TK energy calculation function
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  tfs  --> the time-frequence spectrum (dim: 2D)
 Output: tk   --> the teager energy 
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

    tk(:, ncol) = sqrt(tk_r.^2 + tk_i.^2);
end
tk = max(tk,[],1); % sum(tk, 1);

end

