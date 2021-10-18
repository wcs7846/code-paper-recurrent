function [ fdiv ] = Divergence( evid1, evid2 )
%% DIVERGENCE: this function is used to compute the divergence(belief functions)
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  evid1,evid2 --> the evidence information(a struct array)
 Output: fdiv        --> the divergence between belief function
[Tips]: This function is based on the special design of the system to
calculate some parameters. So this function only can be used under the data
structure in this system.
%}
fdiv = ddiv(evid1, evid2) + ddiv(evid2, evid1);
end

function output = ddiv(evid1, evid2)
%{
 detail: 
 Input:  evid1,evid2 --> the evidence information(a struct array)
 Output: fdiv        --> the divergence between belief function
%}
% qi belongs to evid1; ri belongs to evid2
len1 = length(evid1);
len2 = length(evid2);

output = 0;
for n = 1:len1
    tmp_qi = evid1(n);
    sim_vec = zeros(1, len2); % the similarity vector
    % find the most adapted ri to get the max similarity 
    for m = 1:len2
        tmp_ri = evid2(m);
        if ~isempty(tmp_ri.basicElement) && ~isempty(tmp_qi.basicElement)
            sim_vec(m) = FocalSimilarity(tmp_qi.basicElement, tmp_ri.basicElement);
        else
            sim_vec(m) = 0;
        end
        
    end
    loc = (sim_vec == max(sim_vec));
    
    % ensure to  have only a maximum value
    if sum(loc)>1
        for t = 1:length(loc)
            if loc(t) && t<length(loc)
                loc(t+1:length(loc)) = 0;
            end
        end
    end
    % calculate the ddiv
    % [Ref]:[1]刘娟. 基于信度函数分形维数的模式识别研究[D].西南大学,2014.(Chinese)
    tt = tmp_qi.bpa/(evid2(loc).bpa + eps);
    if tt == 0
        continue;
    end
    output = output + exp(1-max(sim_vec))*tmp_qi.bpa*abs(log(tt));
end
end

function sim = FocalSimilarity(q1, q2)
%{
detail: to compute the similarity between two focal elements
Input:  q1,q2 --> the focal elements
Output: sim   --> the similarity
%}
Attr_List_q1 = extractAttributeList(q1);
Attr_List_q2 = extractAttributeList(q2);
N_attr = length(Attr_List_q1);

mixed_total = 0; union_total = 0;
for n = 1:N_attr
    tmp_q1 = Attr_List_q1{n};
    tmp_q2 = Attr_List_q2{n};
    
    mixed_tmp = calculateMixed(tmp_q1, tmp_q2);
    union_tmp = calculateUnion(tmp_q1, tmp_q2);
    
    mixed_total = length(mixed_tmp) + mixed_total;
    union_total = length(union_tmp) + union_total;
end
sim = mixed_total/union_total;
end

%{
function sim = FocalSimilarity(q1, q2)
%{
detail: to compute the similarity between two focal elements
Input:  q1,q2 --> the focal elements
Output: sim   --> the similarity
%}
Attr_List_q1 = extractAttributeList(q1);
Attr_List_q2 = extractAttributeList(q2);

% compare two attribute list
% <1> num_planes
tmp_q1 = Attr_List_q1.num_planes;
tmp_q2 = Attr_List_q2.num_planes;

mixed_num_planes = calculateMixed(tmp_q1, tmp_q2);
union_num_planes = calculateUnion(tmp_q1, tmp_q2);
% <2> type
tmp_q1 = Attr_List_q1.type;
tmp_q2 = Attr_List_q2.type;

mixed_type = calculateMixed(tmp_q1, tmp_q2);
union_type = calculateUnion(tmp_q1, tmp_q2);
% <3> curvature
tmp_q1 = Attr_List_q1.curvature;
tmp_q2 = Attr_List_q2.curvature;

mixed_curvature = calculateMixed(tmp_q1, tmp_q2);
union_curvature = calculateUnion(tmp_q1, tmp_q2);
% output
mixed_total = length(mixed_num_planes) + length(mixed_type) + length(mixed_curvature);
union_total = length(union_num_planes) + length(union_type) + length(union_curvature);
sim = mixed_total/union_total;
end
%}

function output = calculateMixed(q1, q2)
%{
detail: to calculate the mixed set between two sets.
Input:  q1,q2  --> the atrribute list(cell array)
Output: output --> the mixed set(cell array)
%}
len_q1 = length(q1);
len_q2 = length(q2);
if len_q2 > len_q1
    % swap q1 and q2 to guarantee the q1 longer than q2
    tmp = q1;
    q1 = q2;
    q2 = tmp;
end
% [Tips]: q1 is longer than q2(or equal length)
loc = ismember(q1, q2);

output = q1(loc == 1);
end

function output = calculateUnion(q1, q2)
%{
detail: to calculate the union set between two sets.
Input:  q1     --> the atrribute list(cell array)
Output: output --> the union set(cell array)
%}

tmp = [q1; q2];

% mcs = cellfun(@(x)(mat2str(x)),tmp,'uniformoutput',false);
% [uniqueCells,~,~] = unique(mcs);
% output = uniqueCells;
output = unique(tmp);
end

function output = extractAttributeList(q1)
%{
detail: to extract the attribute list from a focal element
Input:  q1     --> the focal elements
Output: output --> the attribute list(a struct)
%}
N = length(q1);
N_attr = length(q1{1});

output = cell(1,N_attr);
for m = 1:N_attr
    tmp_attr_vec = cell(N,1);
    for n = 1:N
        be = q1{n};
        tmp_attr_vec(n) = be(m);
    end
    attr_vec = unique(tmp_attr_vec);  % used to debug
    output(m) = {attr_vec};
end

end

%{
function output = extractAttributeList(q1)
%{
detail: to extract the attribute list from a focal element
Input:  q1     --> the focal elements
Output: output --> the attribute list(a struct)
%}
[~, N] = size(q1);
output.num_planes = {};
output.type = {};
output.curvature = {};

for n = 1:N
    output.num_planes(n) = {q1(n).num_planes};
    output.type(n) = {abs(q1(n).type)}; % convert char to ascii
    output.curvature(n) = {abs(q1(n).curvature)}; % convert char to ascii
end
% remove the same element 
% <1> num_planes
mcs = cellfun(@(x)(mat2str(x)),output.num_planes,'uniformoutput',false);
[uniqueCells,~,~] = unique(mcs);
output.num_planes = uniqueCells;
% <2> type
mcs = cellfun(@(x)(mat2str(x)),output.type,'uniformoutput',false);
[uniqueCells,~,~] = unique(mcs);
output.type = cellfun(@(x)(char(str2num(x))), uniqueCells, 'uniformoutput',false);
% <3> curvature
mcs = cellfun(@(x)(mat2str(x)),output.curvature,'uniformoutput',false);
[uniqueCells,~,~] = unique(mcs);
output.curvature = cellfun(@(x)(char(str2num(x))), uniqueCells, 'uniformoutput',false);

end
%}
