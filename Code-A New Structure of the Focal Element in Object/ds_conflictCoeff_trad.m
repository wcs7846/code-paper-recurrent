function [ coeff ] = ds_conflictCoeff_trad( evid1, evid2 )
%% DS_CONFLICTCOEFF: this function is used to compute the conflict coefficient(DS evidence)
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  evid1,evid2 --> the evidence information(a struct array)
 Output: coeff       --> the conflict coefficient
[Tips]: This function is based on the special design of the system to
calculate some parameters. So this function only can be used under the data
structure in this system.
%}
len_evid1 = length(evid1);
len_evid2 = length(evid2);
coeff = 0;
% extract the focal elements
for n = 1:len_evid1
    fe1 = evid1(n).basicElement;
    for m = 1:len_evid2
        fe2 = evid2(m).basicElement;
        
        % judge whether it is the same
        if FocalElementMixed(fe1, fe2)
            continue;
        else
            tmp = evid1(n).bpa * evid2(m).bpa;
            coeff = coeff + tmp;
        end
    end
end
coeff = coeff;
end

function res = FocalElementMixed(fe1, fe2)
%{
 detail: to judge the two focal elements have same basic elements
 Input:  fe1,fe2 --> the focal elements(a struct array)
 Output: res     --> true/false
%}
len1 = length(fe1); len2 = length(fe2);
for n = 1:len1
    tmp_be1 = fe1(n);
    for m = 1:len2
        tmp_be2 = fe2(m);
        
        if IsSameBasicElement(tmp_be1, tmp_be2)
            res = true;
            return;
        else
            res = false;
        end
    end
end

end

function res = IsSameBasicElement(be1, be2)
%{
 detail: to judge the two basic elements are same
 Input:  be1,be2 --> the basic elements(a struct array)
 Output: res     --> true/false
%}
N = length(be1);
flag = false(N,1);
for n = 1:N
    if strcmpi(be1{n}, be2{n})
        flag(n) = true;
    end
end
if sum(flag) == N
    res = true;
else
    res = false;
end
end

%{
function res = IsSameBasicElement(be1, be2)
%{
 detail: to judge the two basic elements are same
 Input:  be1,be2 --> the basic elements(a struct array)
 Output: res     --> true/false
%}
flag_num_planes = false;
flag_type = false;
flag_curvatrue = false;

if be1.num_planes == be2.num_planes
    flag_num_planes = true;
end
if be1.type == be2.type
    flag_type = true;
end
if be1.curvature == be2.curvature
    flag_curvatrue = true;
end
res = flag_num_planes & flag_type & flag_curvatrue;
end
%}