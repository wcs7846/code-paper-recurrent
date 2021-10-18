function [ prop ] = pignistic_prop_trans( evid, framework )
%PIGNISTIC_PROP_TRANS 此处显示有关此函数的摘要
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  evid  --> the evidence information(a struct array)
         frame --> the framework
 Output: prop  --> the probability function
[Tips]: This function is based on the special design of the system to
calculate some parameters. So this function only can be used under the data
structure in this system.
%}
N = length(framework.type);e_bpa = 0;
% get the empty set 
for n = 1:length(evid)
    tmp = evid(n).basicElement;
    if isempty(tmp)
        e_bpa = evid(n).bpa; % e_bpa = empty set bpa
        break;
    end
end
% PPT calculation
for n = 1:N
    prop(n).type = framework.type(n);
    prop(n).p = 0;
    % calculate the pignistic probability transformation
    % REF:Smets, P., Kennes, R., 1994. The transferable belief model. Artif. Intell. 66 (2), 191–234.
    N_fe = length(evid);
    for m = 1:N_fe
        fe = evid(m); % select the basic element
        fe_bpa = fe.bpa;
        fe_be  = fe.basicElement;
        N_be = length(fe_be);
        for k = 1:N_be
            if strcmpi(fe_be{k}, prop(n).type)
                prop(n).p = prop(n).p + fe_bpa/length(fe_be) * 1/(1-e_bpa);
                break;
            end
        end
    end
end

end

