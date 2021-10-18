function [ evid_set ] = loadWine_tradition( data, typeClass )
Num_class = length(typeClass);
[N, tt] = size(data);
Num_Attr  = tt-1; % The first is class;
%% create the full library(transform the source data to BPA)
% simple statistics
tt = data(:,1); TotalNum = zeros(Num_class, 1);
for n = 1:Num_class
    TotalNum(n) = length(find(tt == n));
end
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
        if isempty(tmp_vec)
            mean_mat(k,m) = 0;
            std_mat(k,m) = 0;
            continue;
        end
        
        record_mat(k,m) = {tmp_vec};
        % calcualte the mean and std value
        mean_mat(k,m) = mean(tmp_vec);
        std_mat(k,m) = std(tmp_vec);
    end
end

%% generator the BPA for each samples -- traditional structure
evid_set = cell(N, 1);
for k = 1:N % for each sample
    tmp_samples = data(k,2:Num_Attr+1);
    tmp_vec = zeros(Num_Attr, 1);
    fe_vec = cell(Num_class, Num_Attr); % each row is a focal elements
    fe_bpa = zeros(Num_class, Num_Attr);
    
    tt_evid = [];
    for m = 1:Num_Attr
        % generate the BPA for each attributes
        tmp_attr = tmp_samples(m);
        tmp_bpa_vec = zeros(Num_class, 1);
        for n = 1:Num_class
            mu = mean_mat(n,m); s = std_mat(n,m)+eps;
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
        % construct the focal element
        tmp_evid = [];tt = 0; tmp_fe.basicElement = {};
        for n = 1:Num_class
            tmp_fe.basicElement = [tmp_fe.basicElement, fe_vec(n,m)];
            tmp_fe.bpa = fe_bpa(n,m);
            tmp_evid = [tmp_evid, tmp_fe]; tt = tt + tmp_fe.bpa;
        end
        % normalization
        for n = 1:Num_class
            tmp_evid(n).bpa = tmp_evid(n).bpa/tt;
        end
        tt_evid = [tt_evid, {tmp_evid}];
    end
    evid_set(k) = {tt_evid};
end

end

