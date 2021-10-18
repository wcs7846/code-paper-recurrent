function [ cond ] = conductivity( evid_c, inter_dim, aggre_dim)
%% CONDUCTIVITY: a function is used to compute the conductivity of evidence in a evidence set
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  evid_c    --> a evidence set (a cell array)
         inter_dim --> the interaction dimension(evidence and evidence set)
         aggre_dim --> the aggregate dimension(evidence set)
 Output: cond      --> the conductivity

[Tips]: This function is based on the special design of the system to
calculate some parameters. So this function only can be used under the data
structure in this system.
%}

len = length(inter_dim);
loc = inter_dim > aggre_dim;
cond_y1_vec = zeros(1, len);
for n = 1:len
    if loc(n)
        % --> if interaction dimension >  aggregate dimension
        % compute the proj[interaction dim]   
        cap = 1/innerEntropy(evid_c);
        proj_dim = inter_dim(n)*cap;
        
        if proj_dim > aggre_dim
            cond_y1_vec(n) = 0;
        else
            cond_y1_vec(n) = proj_dim;
        end
    else
        % --> if interaction dimension <= aggregate dimension  
        cond_y1_vec(n) = inter_dim(n);
    end    
end
cond = sum(cond_y1_vec);

end

