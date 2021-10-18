function [ innE ] = innerEntropy( evid_c )
%% WITHINENTROPY: this function is used to compute the within entropy of a evidence class
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  evid_c  --> a evidence class (a cell array)
 Output: innE    --> the within entropy of this evidence class
[Tips]: This function is based on the special design of the system to
calculate some parameters. So this function only can be used under the data
structure in this system.
%}
len = length(evid_c);
belief_E_Pyra = zeros(1, len); % belief entropy 
core_E_Pyra   = zeros(1, len); % core entropy
Fdiv_Pyra     = zeros(1, len); % focal divergence

In_E_Pyra     = zeros(1, len); % inner entropy
for n = 1:len
    belief_E_Pyra(n) = beliefEntropy(evid_c{n});
    core_E_Pyra(n)   = coreEntropy(evid_c{n});
    Fdiv_Pyra(n)     = FocalDivergence(evid_c{n});
    
    within_Entropy = belief_E_Pyra(n) + core_E_Pyra(n);
    In_E_Pyra(n) = within_Entropy + Fdiv_Pyra(n);
end

innE = sum(In_E_Pyra);
end

