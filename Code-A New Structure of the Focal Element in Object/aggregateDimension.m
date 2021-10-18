function [ dim ] = aggregateDimension( evid )
%% AGGREGATEDIMENSION: this function is used to compute the aggregate dimension(class)
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  evid  --> a class evidence information(a cell array)
 Output: dim   --> the aggregate dimension of evidence

[Tips]: This function is based on the special design of the system to
calculate some parameters. So this function only can be used under the data
structure in this system.
%}

len = length(evid);
belief_E_Pyra = zeros(1, len); % belief entropy 
core_E_Pyra   = zeros(1, len); % core entropy
Fdiv_Pyra     = zeros(1, len); % focal divergence

In_E_Pyra     = zeros(1, len); % inner entropy
for n = 1:len
    belief_E_Pyra(n) = beliefEntropy(evid{n});
    core_E_Pyra(n)   = coreEntropy(evid{n});
    Fdiv_Pyra(n)     = FocalDivergence(evid{n});
    
    within_Entropy = belief_E_Pyra(n) + core_E_Pyra(n);
    In_E_Pyra(n) = within_Entropy + Fdiv_Pyra(n);
end

Div_mat = zeros(len, len); % divergence matrix of belief functions
for nrow = 1:len
    for ncol = 1:len
        q1 = evid{nrow};
        q2 = evid{ncol};
        
        Div_mat(nrow, ncol) = Divergence(q1, q2);
    end
end
tt = Div_mat - diag(diag(Div_mat)); % remove the diagonal element
Div_Pypr = sum(tt(:))/2 ;
% compute the aggregate dimension of the evidence
dim_Pypr = sum(In_E_Pyra)/(sum(In_E_Pyra)+Div_Pypr); 
dim = dim_Pypr;

end

