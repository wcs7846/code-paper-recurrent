function [ cred_vec, credible_evid ] = credible_degree( Pyra_evid_2 )
%% CREDIBLE_DEGREE: 此处显示有关此函数的摘要
%   此处显示详细说明
len = length(Pyra_evid_2);
%% this way is serial calculation 
% conflict_mat = zeros(len, len);
% for nrow = 1:len
%     for ncol = nrow:len
%         evid1 = Pyra_evid_2{nrow};
%         evid2 = Pyra_evid_2{ncol};
%         tmp = Divergence(evid1, evid2);
% %         tmp= 1 - ds_conflictCoeff(evid1, evid2);
% 
%         conflict_mat(nrow, ncol) = tmp;
%         conflict_mat(ncol, nrow) = tmp;
%     end
% end
%% this way is parallel calculation
t = 1:len;
[nrow,ncol] = meshgrid(t, t);
for t = 1:numel(nrow)
    tmp_conflict_mat(t) = Divergence(Pyra_evid_2{nrow(t)}, Pyra_evid_2{ncol(t)});
end
conflict_mat = reshape(tmp_conflict_mat, len, len);
%% clusting algorithm
Clust = DBSCAN(conflict_mat,mean2(conflict_mat),4);
class_type = unique(Clust);
numOftype = zeros(length(class_type), 1);
for n = 1:length(class_type)
    numOftype(n) = sum(Clust == class_type(n));
end
Trust_class = Clust == class_type(numOftype == max(numOftype));
% mean_conflict = mean(conflict_mat, 2);
% class_mat = zeros(len, len);
% precise = 1e-6;
% for n = 1:len
%     tmp = conflict_mat(n,:);
%     loc = (tmp - mean_conflict(n)) < precise;
%     class_mat(n, loc) = 1;
% end
% class_vec = sum(class_mat, 2);
% Trust_class = class_vec == max(class_vec);
% merge all credible evidence 
% tt = zeros(numel(Trust_class), 1); tt(Trust_class==1) = 1/sum(Trust_class);
% cred_vec = tt;credible_evid = 0;
% return;

credible_evid = 0; flag = 0;
for n = 1:len
    if Trust_class(n) && flag == 0
        credible_evid = Pyra_evid_2{n};
        flag = 1;
        continue;
    end
    
    if Trust_class(n) && flag == 1
        tmp = Pyra_evid_2{n};
        credible_evid = Murphy_evidMerge(credible_evid, tmp);
    end
end
% assume use conflict coeff as credible degree
cred_vec = zeros(len, 1);
for n = 1:len
    cred_vec(n) = 1 - ds_conflictCoeff(credible_evid, Pyra_evid_2{n});
end
cred_vec = cred_vec./sum(cred_vec); 
end

