function [ dim ] = interactionDimension( evid_c, evid )
%% INTERACTIONDIMENSION: this function is used to compute the aggregate dimension(evidence and class)
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  evid_c --> a class evidence information(a cell array)
         evid   --> a unknown class evidence
 Output: dim    --> the interaction dimension between evidence and class

[Tips]: This function is based on the special design of the system to
calculate some parameters. So this function only can be used under the data
structure in this system.
%}

len = length(evid_c);
dim = zeros(1, len);

parfor n = 1:len
    evid_tmp = evid_c{n};
    
    dim(n) = interactDim(evid, evid_tmp);
end

end

function dim = interactDim(evid1, evid2)
%{
detail: to calculate the interaction dimension between two evidence.
Input:  evid1,evid2 --> the evidence 
Output: dim         --> the interaction dimension
%}
belief_e1 = beliefEntropy(evid1);
belief_e2 = beliefEntropy(evid2);

core_e1 = coreEntropy(evid1);
core_e2 = coreEntropy(evid2);

FocalDiv_1 = FocalDivergence(evid1);
FocalDiv_2 = FocalDivergence(evid2);

% within set entropy
within_e1 = belief_e1 + core_e1 + FocalDiv_1;
within_e2 = belief_e2 + core_e2 + FocalDiv_2;

% divergence between two evidence
div_e1e2 = Divergence(evid1, evid2);

dim = (within_e1 + within_e2)/(within_e1 + within_e2 + div_e1e2);
end
