function [ fdiv ] = FocalDivergence( evid )
%% FOCALDIVERGENCE: this function is used to compute the focal divergence(focal elements)
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  evid --> the evidence information(a struct array)
 Output: fdiv --> the focal divergence between focal elements
[Tips]: This function is based on the special design of the system to
calculate some parameters. So this function only can be used under the data
structure in this system.
%}
[~, N] = size(evid);

% handle the bpa == 0 situation in the evidence
tmp_evid = evid; loc = true(N,1);
for n = 1:N
    if tmp_evid(n).bpa == 0
        loc(n) = false;
    end
end
tmp_evid(~loc) = [];
[~,N] = size(tmp_evid);

sim_mat = zeros(N, N);% the similarity matrix 
diff_mat = zeros(N,N);% the difference matrix of BPA
div_mat = zeros(N, N);% the divide matrix of BPA 
for nrow = 1:N
    for ncol = 1:N
        q1 = tmp_evid(nrow).basicElement;
        q2 = tmp_evid(ncol).basicElement;
        
        sim_mat(nrow, ncol)  = FocalSimilarity(q1, q2);
        diff_mat(nrow, ncol) = tmp_evid(nrow).bpa - tmp_evid(ncol).bpa;
        div_mat(nrow, ncol)  = tmp_evid(nrow).bpa/(tmp_evid(ncol).bpa + eps);
    end
end

% calculate the focal divergence 
% [Ref]:[1]刘娟. 基于信度函数分形维数的模式识别研究[D].西南大学,2014.(Chinese)
fdiv_mat = exp(1 - sim_mat).*diff_mat.*log(div_mat);
fdiv = sum(fdiv_mat(:));

end

function sim = FocalSimilarity(q1, q2)
%{
detail: to compute the similarity between two focal elements
Input:  q1,q2 --> the focal elements
Output: sim   --> the similarity
%}
% judge the q1 and q2 is []?
if isempty(q1) && isempty(q2)
    sim = 1;return;
elseif ~isempty(q1) && isempty(q2)
    sim = 0;return;
elseif isempty(q1) && ~isempty(q2)
    sim = 0;return;
end

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
%{
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
%}
%{
function output = calculateUnion(q1, q2)
%{
detail: to calculate the union set between two sets.
Input:  q1     --> the atrribute list(cell array)
Output: output --> the union set(cell array)
%}

tmp = [q1, q2];

mcs = cellfun(@(x)(mat2str(x)),tmp,'uniformoutput',false);
[uniqueCells,~,~] = unique(mcs);
output = uniqueCells;
end
%}