% 这个脚本是用来显示不同葡萄糖浓度下的光声信号
clear; close all;

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
mu_a_data = [];
absorption_data = [];
for n = 1:len
    if isempty(newStr)
        continue;
    end
    fileName = strcat(head_ac, num2str(gc(n)), '.mat');
    path_tmp = strcat(path, '\', fileName);
    
    load(path_tmp); % load the data
    
    mu_a = vmcmedium.absorption_coefficient;
    [row,col] = size(mu_a);
    sample_index = round(row/2);
    data = mu_a(sample_index, :)';
    mu_a_data = [mu_a_data, data];
    
%     sensor_data = sensor_data';
%     [row, num_probe] = size(sensor_data);
%     data = sensor_data(:, round(num_probe/2));
    
    [row_li, col] = size(solution.grid_fluence);
    data = vmcmedium.absorbed_energy(round(row_li/2), :)';
    absorption_data = [absorption_data, data];
end

%% show the data
mu_a_show = figure;
light_show = figure;
dt = 0.01*1e3; 
xrange = (1:row)*dt; lineColor = linspecer(len); 
label = [];
for n = 1:len  
    show_x = find(xrange > 3800 & xrange < 4200 == 1) ;
    
    label = [label, string(strcat(num2str(gc(n)), 'mg/L'))]; % string function must be used on MATLAB2016b
    
    show_data = mu_a_data(:, n);
    figure(mu_a_show);plot(xrange(show_x), show_data(show_x), 'color', lineColor(n, :)); hold on;
    show_data = absorption_data(:, n);
    figure(light_show);plot(xrange(show_x), show_data(show_x), 'color', lineColor(n, :)); hold on;
end

figure(mu_a_show);
colorbar('Ticks',0:1/(length(gc)-1):1,'TickLabels', label);colormap(lineColor);
xlabel('depth/um'); ylabel('mu_a / mm^-^1'); title(strcat('the absorption coefficient distribution', ' (fs:', num2str(round(1/dt)) ,'Hz)'));
xlim([min(xrange(show_x)), max(xrange(show_x))]);

figure(light_show);
colorbar('Ticks',0:1/(length(gc)-1):1,'TickLabels', label);colormap(lineColor);
xlabel('depth/um'); ylabel('absorption energy'); title(strcat('absorption energy data', ' (fs:', num2str(round(1/dt)) ,'Hz)'));
xlim([min(xrange(show_x)), max(xrange(show_x))]);


