function [ bpa, save_interDim ] = BPA_Project( evid_c, evid_un, aggre_dim )
%% BPA_PROJECT
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  evid_c    --> all evidence sets (a cell array for all evidence set)
         evid_un   --> the unknown class evidence
         aggre_dim --> the aggregate dimension(for all evidence set)
 Output: bpa       --> the BPA based on conductivity
         save_interDim --> used to save the interDim
[Tips]: This function is based on the special design of the system to
calculate some parameters. So this function only can be used under the data
structure in this system.
%}
% evid_c = [Pyramid, Lshape, Handle, Cylinder];
Pyra_evid = evid_c{1};
Lshape_evid = evid_c{2};
Handle_evid = evid_c{3};
Cylin_evid = evid_c{4};

InterDim_Pyra = interactionDimension(Pyra_evid, evid_un);     % Pyramid evidence
InterDim_Lshape = interactionDimension(Lshape_evid, evid_un); % L shape evidence
InterDim_Handle = interactionDimension(Handle_evid, evid_un); % Handle evidence
InterDim_Cylin = interactionDimension(Cylin_evid, evid_un);   % Cylinder evidence
% save the intermediate result
save_interDim.Pyra   = InterDim_Pyra;
save_interDim.Lshape = InterDim_Lshape;
save_interDim.Handle = InterDim_Handle;
save_interDim.Cylin  = InterDim_Cylin;
% compute the conductivity of Y1 in class Pyramid
% aggre_dim = [Pyramid, Lshape, Handle, Cylinder];
weight = 1- (abs(InterDim_Pyra - aggre_dim(1))/aggre_dim(1));
proj_dim_Pyra = weight.*InterDim_Pyra;
cond_y1_Pyra = mean(proj_dim_Pyra);

weight = 1- (abs(InterDim_Lshape - aggre_dim(2))/aggre_dim(2));
proj_dim_Lshape = weight.*InterDim_Lshape;
cond_y1_Lshape = mean(proj_dim_Lshape);

weight = 1- (abs(InterDim_Handle - aggre_dim(3))/aggre_dim(3));
proj_dim_Handle = weight.*InterDim_Handle;
cond_y1_Handle = mean(proj_dim_Handle);

weight = 1- (abs(InterDim_Cylin - aggre_dim(4))/aggre_dim(4));
proj_dim_Cylin = weight.*InterDim_Cylin;
cond_y1_Cylin = mean(proj_dim_Cylin);

y1_BPA = [cond_y1_Pyra, cond_y1_Lshape, cond_y1_Handle, cond_y1_Cylin];
bpa = y1_BPA/sum(y1_BPA);

end

