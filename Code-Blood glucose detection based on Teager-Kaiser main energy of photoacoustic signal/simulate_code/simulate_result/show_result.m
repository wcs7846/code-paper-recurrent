% 这个脚本是用来显示不同葡萄糖浓度下的光声信号
clear;close all;

path = '.\acoustic';
list = dir(path);
list = list(3:length(list));
% check
if isempty(list)
    error("The acoustic folder is empty!");
end

% extract the concentration of glucose
len = length(list);
gc = zeros(len, 1);

head_ac = 'acoustic_skin_'; % the acoustic file header
for n = 1:len
    % sperate the name 
    folder_name = list(n).name;    
    newStr = folder_name(length(head_ac)+1:length(folder_name)-4);
    %
    gc(n) = eval(newStr);
end
gc = sort(gc); % sort the number

%% load the data
% Tips: only extract the center probe data for every concentration
acoustic_data = [];
absorption_data = [];
for n = 1:len
    if isempty(newStr)
        continue;
    end
    fileName = strcat(head_ac, num2str(gc(n)), '.mat');
    path_tmp = strcat(path, '\', fileName);
    
    load(path_tmp); % load the data
    
    sensor_data = sensor_data';
    [row_ac, num_probe] = size(sensor_data);
    data = sensor_data(:, round(num_probe/2));
    acoustic_data = [acoustic_data, data];
    
    [row_li, col] = size(solution.grid_fluence);
    data = vmcmedium.absorbed_energy(round(row_li/2), :)';
    absorption_data = [absorption_data, data];
end

%% show the data
acoustic_show = figure; light_show = figure;
dt = 2.5e-9; % this data can see the skin_model_simulation.m(kgrid.dt)
xrange = (1:row_ac); lineColor = linspecer(len); 
label = [];
for n = 1:len  
    label = [label, string(strcat(num2str(gc(n)), 'mg/L'))]; % string function must be used on MATLAB2016b
    
    figure(acoustic_show);plot(xrange*dt, acoustic_data(:, n), 'color', lineColor(n, :)); hold on;
    figure(light_show);plot((1:row_li)*0.01, absorption_data(:, n), 'color', lineColor(n, :)); hold on;
end
figure(acoustic_show);legend(label);
xlabel('time/s'); ylabel('pressure'); title(strcat('the center porbe acoustic data', ' (fs:', num2str(round(1/dt)) ,'Hz)'));

figure(light_show);legend(label);
xlabel('depth/mm'); ylabel('absorption energy'); title(strcat('absorption energy data', ' (fs:', num2str(round(1/dt)) ,'Hz)'));

sensor_data = acoustic_data'; 
save('acoustic_skin','sensor_data');