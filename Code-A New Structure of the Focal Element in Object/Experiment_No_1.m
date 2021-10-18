%% Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
% This script is designed to test the information fractals for evidential
% pattern classification
% Author: MarkLHF
% Date: 2020-3-13(first version)
clear; clc; close all;
%% create the full library
lib.num_planes_select = 1:8;
% [Tips]: must convert the cell array to char array
lib.type_select = char({'t', 's', 'r', 'b', 'c', 'd'}); 
lib.curvature_select = char({'p', 'q'});

total_num = length(lib.num_planes_select) * ...
            length(lib.curvature_select) * ...
            length(lib.type_select);
        
Know_lib = createKnowledgeLib(lib);
%% class: Pyramid
% vector: {(4,T,P) (1,S,P) (1,R,P) FF}
p_4tp = createBasicElement(4, 't', 'p');
p_1sp = createBasicElement(1, 's', 'p');
p_1rp = createBasicElement(1, 'r', 'p');

pf_4tp = createFocalElement(p_4tp, 0.5);
pf_1sp = createFocalElement(p_1sp, 0.3);
% --> generate a1
a1 = createEvidence([pf_4tp, pf_1sp], Know_lib);

pf_1rp = createFocalElement(p_1rp, 0.25);
% --> generate a2
a2 = createEvidence([pf_4tp, pf_1rp], Know_lib);

pf_4tp = createFocalElement(p_4tp, 0.7);
% --> generate a3
a3 = createEvidence(pf_4tp, Know_lib);

% p_2cp = createBasicElement(2, 'c', 'p');
% p_1rc = createBasicElement(1, 'r', 'c');
% pf_2cp = createFocalElement(p_2cp, 0.4);
% pf_1rc = createFocalElement(p_1rc, 0.4);
% a4 = createEvidence([pf_2cp, pf_1rc], Know_lib);

p_6rp = createBasicElement(6, 'r', 'p');
pf_6rp = createFocalElement(p_6rp, 0.5);
a4 = createEvidence(pf_6rp, Know_lib);

p_2bp = createBasicElement(2, 'b', 'p');
p_2sp = createBasicElement(2, 's', 'p');
pf_2bp = createFocalElement(p_2bp, 0.4);
pf_6rp = createFocalElement(p_6rp, 0.2);
pf_2sp = createFocalElement(p_2sp, 0.2);
a5 = createEvidence([pf_2bp, pf_6rp, pf_2sp], Know_lib);

% create a evidence cell array
% <1> Pyramid
Pyra_evid_1 = {a1, a2, a3};
Pyra_evid_2 = {a1, a2, a3, a4, a5};
%% generate 3 types evidence
p_2tp = createBasicElement(2, 't', 'p');
p_5rp = createBasicElement(5, 'r', 'p');
p_6rp = createBasicElement(6, 'r', 'p');

pf_4tp = createFocalElement(p_4tp, 0.4);
pf_1sp = createFocalElement(p_1sp, 0.4);
evid_r = createEvidence([pf_4tp, pf_1sp], Know_lib); % right evidence

pf_2tp = createFocalElement(p_2tp, 0.4);
pf_1sp = createFocalElement(p_1sp, 0.2);
pf_5rp = createFocalElement(p_5rp, 0.1);
evid_pr = createEvidence([pf_2tp, pf_1sp, pf_5rp], Know_lib); % partly right evidence

pf_6rp = createFocalElement(p_6rp, 0.5);
evid_w = createEvidence(pf_6rp, Know_lib); % wrong evidence
%% compute the correlation 

% aggregateDimension(Pyra_evid_1);
% aggregateDimension(Pyra_evid_2);

cred_vec = credible_degree(Pyra_evid_2);


