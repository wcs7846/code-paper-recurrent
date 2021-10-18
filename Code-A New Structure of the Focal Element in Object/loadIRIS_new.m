function [ evid_set ] = loadIRIS_new( data, TotalNum, typeClass )
%% basic information
Num_class = length(typeClass);
[N, tt] = size(data);
Num_Attr  = tt-1; % The first is class;
%% create the full library(transform the source data to BPA)
% calculate the mean and std of every attributes of every class
mean_mat = zeros(Num_class, Num_Attr); % row: same class; col: same attribute
std_mat  = zeros(Num_class, Num_Attr);
record_mat = cell(Num_class, Num_Attr);
for m = 1:Num_Attr % for each attribute
    
    for k = 1:Num_class
        tmp_vec = zeros(TotalNum(k),1); t=1;
        for n = 1:N % for each sample
            tmp_attr = data(n,m+1);
            if k == data(n,1)
                tmp_vec(t) = tmp_attr;
                t = t+1;
            end
        end
        record_mat(k,m) = {tmp_vec};
        % calcualte the mean and std value
        mean_mat(k,m) = mean(tmp_vec);
        std_mat(k,m) = std(tmp_vec);
    end
end

%% generator the BPA for each samples -- new structure
evid_set = cell(N, 1);
for k = 1:N % for each sample
    tmp_samples = data(k,2:Num_Attr+1);
    tmp_vec = zeros(Num_Attr, 1);
    fe_vec = cell(Num_class, Num_Attr); % each row is a focal elements
    fe_bpa = zeros(Num_class, Num_Attr);
    for m = 1:Num_Attr
        % generate the BPA for each attributes
        tmp_attr = tmp_samples(m);
        tmp_bpa_vec = zeros(Num_class, 1);
        for n = 1:Num_class
            mu = mean_mat(n,m); s = std_mat(n,m);
            tmp_func = @(x)(1/sqrt(2*pi*s^2)*exp(-(x-mu)^2/(2*s^2)));
            
            tt = tmp_func(tmp_attr);
            if (tmp_attr<mu-3*s)||(tmp_attr>mu+3*s)
                tt = 0;
            end
            tmp_bpa_vec(n) = tt;
        end
        tmp_bpa_vec = tmp_bpa_vec/(sum(tmp_bpa_vec)+eps); % for No.44
        tt = [tmp_bpa_vec, (1:Num_class)'];
        
        % allocate the basic element
        tt2 = sortrows(tt, 1, 'descend');
        for n = 1:Num_class
            fe_vec(n,m) = typeClass(tt2(n,2));
            fe_bpa(n,m) = tt2(n,1);
        end
    end
    % construct the focal element
    tmp_evid = [];tt = 0;uni_set = createKnowledgeLib(typeClass, Num_Attr);
    for n = 1:Num_class
        tmp_bpa = mean(fe_bpa(n,:));
        tmp_fe.basicElement = {fe_vec(n,:)};
        tmp_fe.bpa = tmp_bpa;
        tmp_evid = [tmp_evid, tmp_fe]; tt = tt + tmp_bpa;
    end
    % for universal set
    tmp_fe.basicElement = uni_set;
    if tt > 1
        for n = 1:Num_class
            tmp_evid(n).bpa = tmp_evid(n).bpa/tt;
        end
        tmp_fe.bpa = 0;
    else
        tmp_fe.bpa = 1-tt;
    end
    tmp_evid = [tmp_evid, tmp_fe];
    evid_set(k) = {tmp_evid};
end

end

