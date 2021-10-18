function [ merge_evid ] = Murphy_evidMerge_trad( evid1, evid2 )
%% DS_EVIDMERGE: this function is used to compute the evidence merge(DS evidence)
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  evid1,evid2 --> the evidence information(a struct array)
 Output: merge_evid  --> the merge evidence
[Tips]: This function is based on the special design of the system to
calculate some parameters. So this function only can be used under the data
structure in this system.
%}

% compute all focal elements in the merge evidence
merge_fe = mergeFocalElements(evid1, evid2);
len = length(merge_fe);

% murphy average rules
for n = 1:len
    len1 = length(evid1); len2 = length(evid2);
    tmp_fe = merge_fe(n).basicElement;
    tmp_evid(n).basicElement = merge_fe(n).basicElement;
    
    bpa_evid1 = 0;
    for m = 1:len1
        tmp_evid1_fe = evid1(m).basicElement;
        if IsSameFE(tmp_fe, tmp_evid1_fe)
            bpa_evid1 = evid1(m).bpa;
            break;
        end
    end
    
    bpa_evid2 = 0;
    for m = 1:len2
        tmp_evid2_fe = evid2(m).basicElement;
        if IsSameFE(tmp_fe, tmp_evid2_fe)
            bpa_evid2 = evid2(m).bpa;
            break;
        end
    end
    
    tmp_evid(n).bpa = 1/2*(bpa_evid1 + bpa_evid2);
end

% compute the conflict coefficient
K = ds_conflictCoeff_trad(tmp_evid, tmp_evid);
% compute the coressponding BPA 
for n = 1:len
    merge_fe(n).bpa = 1/(1-K+eps)*mergeBPA(tmp_evid, tmp_evid, merge_fe(n).basicElement);
end

merge_evid = merge_fe;
end

function merge_bpa = mergeBPA(evid1, evid2, fe)
%{
 detail: to compute all focal elements in the merge evidence
 Input:  evid1,evid2 --> the evidences(a struct array)
         fe          --> the focal element
 Output: merge_bpa   --> the corresponding bpa
%}
len1 = length(evid1); len2 = length(evid2);
merge_bpa = 0;
for n = 1:len1
    fe1 = evid1(n).basicElement;
    for m = 1:len2
        fe2 = evid2(m).basicElement;
        
        InterSec = intersectionFocalElements(fe1, fe2);
        
        if IsSameFE(InterSec, fe)
            merge_bpa = merge_bpa + evid1(n).bpa * evid2(m).bpa;
        end
    end
end

end

function fe = intersectionFocalElements(fe1, fe2)
%{
 detail: to compute the intersection part between fe1 and fe2
 Input:  fe1,fe2 --> the evidences(a struct array)
 Output: fe      --> the intersection part
%}
fe = [];
len1 = length(fe1); len2 = length(fe2);
% make sure that the tmp_fe1 is not less than tmp_fe2(tmp_fe1 > tmp_fe2)
if len2 > len1
    tmp_fe1 = fe2;
    tmp_fe2 = fe1;
else
    tmp_fe1 = fe1;
    tmp_fe2 = fe2;
end
% The criterion of computing the same part of two focus elements
% --> 1. extract all basic elements of fe2 (named as be2)
% --> 2. search whether the same focal element exists for fe2
for m = 1:length(tmp_fe2)
    be2 = tmp_fe2(m);
    flag_vec = zeros(1, length(tmp_fe1));
    for n = 1:length(tmp_fe1)
        be1 = tmp_fe1(n);
        if IsSameBasicElement(be1, be2)
            flag_vec(n) = true;
        end
    end
     
    if sum(flag_vec) > 0
        % add the focal elements into the output
        fe =[fe, tmp_fe2(m)];
    end  
end

end

function merge_fe = mergeFocalElements(evid1, evid2)
%{
 detail: to compute all focal elements in the merge evidence
 Input:  evid1,evid2 --> the evidences(a struct array)
 Output: merge_fe    --> true/false
%}
len1 = length(evid1); len2 = length(evid2);
%
for n = 1:len1
    merge_fe(n).basicElement = evid1(n).basicElement;
end
% find the focal elements in the evid2 and not in the evid1
tt = 1;
for m = 1:len2
    fe2 = evid2(m).basicElement;
    flag_vec = zeros(1, len1); % to record the states between fe2 and evid1
    for n = 1:len1
        fe1 = evid1(n).basicElement;
        
        if IsSameFE(fe1, fe2)
            flag_vec(n) = 1;            
        end
    end
    % judge whether the fe2 need to be added in merge evidence
    if sum(flag_vec) == 0
        % add the focal elements into the merge evidence
        merge_fe(len1+tt).basicElement = evid2(m).basicElement;
        tt = tt + 1;
    end  
end
end

function res = IsSameFE(fe1, fe2)
%{
 detail: to judge whether two focal elements are same
 Input:  fe1,fe2 --> the focal elements(a struct array)
 Output: res     --> true/false
%}
len1 = length(fe1); len2 = length(fe2);
if len1 ~= len2 
    res = false; return;
end
% make sure that the tmp_fe1 is not less than tmp_fe2
if len2 > len1
    tmp_fe1 = fe2;
    tmp_fe2 = fe1;
else
    tmp_fe1 = fe1;
    tmp_fe2 = fe2;
end
% The criterion of judging the same two focus elements
% --> 1. extract all basic elements of fe2 (named as be2)
% --> 2. search whether the same focal element exists for fe2
for m = 1:length(tmp_fe2)
    be2 = tmp_fe2(m);
    flag_vec = zeros(1, length(tmp_fe1));
    for n = 1:length(tmp_fe1)
        be1 = tmp_fe1(n);
        if IsSameBasicElement(be1, be2)
            flag_vec(n) = true;
        end
    end
    
    if sum(flag_vec)
        res = true;
    else
        res = false;
        return;
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