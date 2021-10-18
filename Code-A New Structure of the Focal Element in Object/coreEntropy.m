function [ corE ] = coreEntropy( evid )
%% COREENTROPY: this function is used to compute the core entropy of a evidence
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  evid --> the evidence information(a struct array)
 Output: corE --> the core entropy of this evidence
[Tips]: This function is based on the special design of the system to
calculate some parameters. So this function only can be used under the data
structure in this system.
%}
[~,N] = size(evid);

tmp_evid = evid; loc = true(N,1);
for n = 1:N
    if tmp_evid(n).bpa == 0
        loc(n) = false;
    end
end
tmp_evid(~loc) = [];

[~,N] = size(tmp_evid);
h_vec = zeros(1, N);

for n = 1:N
    [~, h_vec(n)] = size(tmp_evid(n).basicElement); 
end
h_vec(h_vec == 0) = []; % remove the [] focal elements


total_FE = sum(h_vec);
h_vec = h_vec/total_FE;

corE = -sum(h_vec.*log(h_vec));
end

