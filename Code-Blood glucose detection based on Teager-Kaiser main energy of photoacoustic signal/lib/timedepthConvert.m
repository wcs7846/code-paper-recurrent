function [ depth_vec, signal_valid ] = timedepthConvert( signal_vec, medium_speed, medium_thick, dt)
%% TIMEDEPTHCONVERT: time depth conversion
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  signal_vec   --> the time series of signal
         medium_speed --> the sound speed in the medium (unit: m/s)
         medium_thick --> the thickness of the medium   (unit: um)
         dt           --> the time interval
 Output: depth_vec    --> the depth series (unit: um)
         signal_valid --> the clipped signal (remove the invalid part)
%}
time_used = medium_thick ./ (medium_speed*1e6);
step_used = round(time_used / dt);

depth_vec = zeros(1, sum(step_used));
t_left = 1; t_right = 0; tt = 0;
for n = 1:length(step_used)
    t_left = 1 + t_right;
    t_right = step_used(n) + t_right;
    
    depth_vec(t_left:t_right) = (1:step_used(n))*dt*medium_speed(n)*1e6 + tt;
    tt = depth_vec(t_right);
end
signal_valid = signal_vec(1:sum(step_used));
end

