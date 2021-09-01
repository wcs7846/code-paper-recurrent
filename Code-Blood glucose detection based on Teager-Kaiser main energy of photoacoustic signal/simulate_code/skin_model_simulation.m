%% Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
% This script is designed to simulate the photoacoustic sensor recoder
% Author: MarkLHF
% Date: 2019-9-27

% ToolBox: ValoMC + k-Wave
% clear;; clc;
close all;
%% [1]Target: an optical field simulation of two-dimensional layer media
% setting the toolbox
addpath(genpath('D:\Matlab_program\Matlab_toolbox\ValoMC')); % ValoMC toolbox
addpath(genpath('D:\Matlab_program\Matlab_toolbox\k-Wave')); % k-Wave toolbox
disp("Starting to simulate the optical transmission");
%% 1. create a ValoMC mesh
% the skin model is 3-layer medium and 1D medium (lamda: 1064nm)
layer_e = 0.1; %[mm] epidermis
layer_d = 1.4; %[mm] dermis
layer_s = 2.0; %[mm] subcutaneous fat
layer_m = 2.0; %[mm] muscle
% y axis
y = (layer_e+layer_d+layer_s+layer_m);  %[mm]
dh = 0.01;% grid point spacing in the x/y direction [mm] = 10 um
% Create the k-Wave grid
dx = dh/1e3;        % grid point spacing in the x direction [m]
dy = dh/1e3;        % grid point spacing in the y direction [m]
Ny = round((layer_e+layer_d+layer_s+layer_m)/dh);     % number of grid points in the x (row) direction
Nx = round(y/dh);           % number of grid points in the y (column) direction
kgrid = makeGrid(Nx, dx, Ny, dy);
vmcmesh = createGridMesh(kgrid.y_vec*1e3, kgrid.x_vec*1e3); % [m to mm]
%% 2. Define optical coefficients
% [2.1] initialization background coefficients
vmcmedium.absorption_coefficient = zeros(Nx,Ny); % absorption coefficient [1/mm]
vmcmedium.scattering_coefficient = zeros(Nx,Ny); % scattering coefficient [1/mm]
vmcmedium.scattering_anisotropy  = zeros(Nx,Ny); % scattering anisotropy parameter [unitless]
vmcmedium.refractive_index = zeros(Nx,Ny);       % refractive index [unitless]
%% epidermis layer proprety setting(as the background)
depth_range = 1:round(layer_e/dh);
mua_e = 0.04; mus_e = 3; % the absorption and scattering coefficient [1/mm]
g_e = 0.8 ; n_e = 1.335;

vmcmedium.absorption_coefficient(:, depth_range) = mua_e;
vmcmedium.scattering_coefficient(:, depth_range) = mus_e;
vmcmedium.scattering_anisotropy(:, depth_range)  = g_e;
vmcmedium.refractive_index(:, depth_range)       = n_e;
%% dermis layer proprety setting
depth_range = round(layer_e/dh)+1 : round((layer_e+layer_d)/dh);
mua_d = 0.06; mus_d = 1.8; % the absorption and scattering coefficient [1/mm]
g_d = 0.8 ; n_d = 1.37;

vmcmedium.absorption_coefficient(:, depth_range) = mua_d;
vmcmedium.scattering_coefficient(:, depth_range) = mus_d;
vmcmedium.scattering_anisotropy(:, depth_range)  = g_d;
vmcmedium.refractive_index(:, depth_range)       = n_d;
%% subcutanuous fat layer proprety setting
depth_range = round((layer_e+layer_d)/dh)+1 : round((layer_e+layer_d+layer_s)/dh);
mua_s = 0.1; mus_s = 1.6; % the absorption and scattering coefficient [1/mm]
g_s = 0.8 ; n_s = 1.4;

vmcmedium.absorption_coefficient(:, depth_range) = mua_s;
vmcmedium.scattering_coefficient(:, depth_range) = mus_s;
vmcmedium.scattering_anisotropy(:, depth_range)  = g_s;
vmcmedium.refractive_index(:, depth_range)       = n_s;
%% muscle layer proprety setting
depth_range = round((layer_e+layer_d+layer_s)/dh)+1 : round((layer_e+layer_d+layer_s+layer_m)/dh);
mua_m = 0.035; mus_m = 0.6; % the absorption and scattering coefficient [1/mm]
g_m = 0.9 ; n_m = 1.4;

vmcmedium.absorption_coefficient(:, depth_range) = mua_m;
vmcmedium.scattering_coefficient(:, depth_range) = mus_m;
vmcmedium.scattering_anisotropy(:, depth_range)  = g_m;
vmcmedium.refractive_index(:, depth_range)       = n_m;
%% set blood and blood vessal proprety 
mua_blood = blood_model(0.5); % assume 1 g/L (normal: 0.07~0.13g/L)
mus_blood = 64.2;
g_blood = 0.98; n_blood = 1.36;
depth_loc = 4.0; % [mm]
width = 0.1;     % [mm] = 0.1mm to 10mm  

width_n = round(width/dh);
depth_range = round(depth_loc/dh) - width_n + 1 : round(depth_loc/dh) + width_n;

vmcmedium.absorption_coefficient(:, depth_range) = mua_blood;
vmcmedium.scattering_coefficient(:, depth_range) = mus_blood;
vmcmedium.scattering_anisotropy(:, depth_range)  = g_blood;
vmcmedium.refractive_index(:, depth_range)       = n_blood;
%% 3. Define light source
% Set a light source 
vmcboundary = createBoundary(vmcmesh);
line_start = [-3*(layer_e+layer_d+layer_s+layer_m)/4, 0];
line_end   = [0,0];
line_width = 50*dh;

lightsource = findBoundaries(vmcmesh, 'direction', ...
                                           line_start, ...
                                           line_end, ...
                                           line_width); 
% set the location of light source      
% vmcboundary.lightsource(lightsource) = {'direct'};
vmcboundary.lightsource(lightsource) = {'gaussian'}; % cosinic
vmcboundary.lightsource_gaussian_sigma(lightsource) = 0.01;
% set the direction of light source --> setting by yourself
lightsource_direction = line_end - line_start;
vmcboundary.lightsource_direction(lightsource, 1) = lightsource_direction(1);
vmcboundary.lightsource_direction(lightsource, 2) = lightsource_direction(2);
vmcboundary.lightsource_direction_type(lightsource) = {'absolute'}; 
%% Run the Monte Carlo simulation
options.photon_count = 2e7; % default is 1e6
solution = ValoMC(vmcmesh, vmcmedium, vmcboundary, options);

%% show 
figure;
handle = patch('Faces',vmcmesh.H,'Vertices',vmcmesh.r,'FaceVertexCData', solution.element_fluence, 'FaceColor', 'flat','EdgeColor','none');
xlabel('y-position [mm]');
ylabel('x-position [mm]');
colormap default;

c = colorbar;  % create a colorbar
axis image;
title('grid fluence [J/m3]');
% plot the fluence near the source position
% figure;
% y_range = [-5, 5]; % [mm]
% y_range_index = y_range(1)/dy : y_range(2)/dy;
% plot(y_range_index*dy, solution.grid_fluence(1, y_range_index+round(Ny/2)), 'b-');hold on;
%% [2]Target: calculate the absorption energy distribution

% Define the Gruneisen parameter describing photoacoustic efficiency
vmcmedium.gruneisen_parameter = 0.12*ones(Nx, Ny);      % [unitless]
%% Compute the initial pressure from the photon fluence
% Note that since the medium was defined as two dimensional arrays,
% the output is also given as a two-dimensional array.

% Compute the absorbed optical energy density
% 1e3 converts [J/mm^2] to [J/m^2]
vmcmedium.absorbed_energy = vmcmedium.absorption_coefficient .* solution.grid_fluence*1e3; % [J/m3]
%% show the absorbed energy distribution
figure('rend','painters');
imagesc(kgrid.y_vec*1e3, kgrid.x_vec*1e3, vmcmedium.absorbed_energy);
xlabel('y-position [mm]');
ylabel('x-position [mm]');
colormap default; colorbar;
axis image;
title('absorbed energy [J/m3]');

source.p0 = vmcmedium.gruneisen_parameter .* vmcmedium.absorbed_energy;  % [Pa]
absorb_energy    = vmcmedium.absorbed_energy;

save_path = './simulate_result/optical/';
save(strcat(save_path, 'optical_skin.mat'),'kgrid', 'source', 'absorb_energy', 'vmcmedium', 'solution');
% return;
%% [3]Target: a ultrasonic field simulation of two-dimensional layer media
disp("Starting to simulate the ultrasonic transmission");
medium.sound_speed = zeros(Nx, Ny);    % [m/s]
medium.density     = zeros(Nx, Ny);    % [kg/m^3]
%% epidermis layer proprety setting(as the background)
depth_range = 1:round(layer_e/dh);
den_e = 1200; c_e = 1645; % the density [kg/m^3] and sound speed [m/s]

medium.sound_speed(:, depth_range) = c_e;
medium.density(:, depth_range)     = den_e;
%% dermis layer proprety setting
depth_range = round(layer_e/dh)+1 : round((layer_e+layer_d)/dh);
den_d = 1200; c_d = 1595; % the density and sound speed

medium.sound_speed(:, depth_range) = c_d;
medium.density(:, depth_range)     = den_d;
%% subcutanuous fat layer proprety setting
depth_range = round((layer_e+layer_d)/dh)+1 : round((layer_e+layer_d+layer_s)/dh);
den_s = 1000; c_s = 1450; % the density and sound speed

medium.sound_speed(:, depth_range) = c_s;
medium.density(:, depth_range)     = den_s;
%% muscle layer proprety setting
depth_range = round((layer_e+layer_d+layer_s)/dh)+1 : round((layer_e+layer_d+layer_s+layer_m)/dh);
den_m = 1050; c_m = 1547; % the density and sound speed

medium.sound_speed(:, depth_range) = c_m;
medium.density(:, depth_range)     = den_m;
%% set blood and blood vessal proprety 
den_blood = 1000; c_blood = 1500; % the density and sound speed
depth_loc = 4.0; % [mm]
width = 0.1;     % [mm] = 0.1mm to 10mm  

width_n = round(width/dh);
depth_range = round(depth_loc/dh) - width_n + 1 : round(depth_loc/dh) + width_n;
medium.sound_speed(:, depth_range) = c_blood;
medium.density(:, depth_range)     = den_blood;
%% 
% l_tmp = round(L1/dh);
% medium.sound_speed(disc_indices) = 1800;       % [m/s]
% medium.density(disc_indices) = 1400;       % [m/s]
% Define a linear sensor
width = 20; % even number
startpoint = [round(Ny/2-width/2), 1];  % [grid points]
endpoint   = [round(Ny/2+width/2), 1];  % [grid points]

sensor.mask = makeLine(Nx, Ny, startpoint, endpoint);
%% Move the perfectly matched layer (PML) outside of the computation domain and run the acoustic simulation
% PML outside       --> 'PMLInside' == false
% PML invisibility  --> 'PlotPML'   == false
% PML absorption    --> 'PMLAlpha' default 2
fs = 5e7; % 50MHz
% kgrid.dt = 1/fs;
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false, 'PlotPML', false, 'PMLAlpha', 5);
%% show 
% plot the simulated sensor data
figure;
mesh(sensor_data);
title('Sensor data');
hold off
% initial_pressure = source.p0;
% 
figure;
len = length(sensor_data(round(width/2),:));
plot(kgrid.dt*(1:len), sensor_data(round(width/2),:), 'b-.');hold on;
xlabel('time/s');ylabel('pressure');
% sensor_data = sensor_data';
save_path = './simulate_result/acoustic/';
save(strcat(save_path, 'acoustic_skin.mat'), 'source', 'absorb_energy', 'sensor_data', 'vmcmedium', 'solution', 'kgrid');
% save('simulate_full.mat', 'initial_pressure', 'absorb_energy', 'sensor_data');