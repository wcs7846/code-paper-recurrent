function [ belE ] = beliefEntropy( evidence )
%% BELIEFENTROPY: this function is used to compute the belief entropy of a evidence
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  evid --> the evidence information(a struct array)
 Output: belE --> the belief entropy of this evidence
[Tips]: This function is based on the special design of the system to
calculate some parameters. So this function only can be used under the data
structure in this system.
%}
[~, N] = size(evidence);

belE = 0;
for n = 1:N
    tmp = evidence(n).bpa;
    
    if tmp == 0
        continue;
    end
    belE = belE + tmp*log(tmp);
end
belE = -belE;
end

