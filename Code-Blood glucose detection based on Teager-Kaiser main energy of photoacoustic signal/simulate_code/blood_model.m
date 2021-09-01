function [ output_args ] = blood_model( c )
%% BLOOD_MODEL calculate the absorption coefficient for different
% concentration glucose solution
%{
 Copyright(C),UESTC, School of Information and Communication Engineering, IDIP
 detail
 Input:  c_glucose  --> the concentration (unit: g/L)
 Output: output_args --> the absorption coefficient(unit: mm-1)
%}

%% assume the four main component in the blood 
% blood = water + HGb(hemoglobin) + glucose + lipid
% the Combination Rule is ua = ln(10)(ci*si)
% ci : the molar concentration (unit: mol/L)
% si : the molar extinction coefficient (unit: L*mol-1*cm-1)

%% the molar concentration (unit: mol/L)
c_HGb = 2.32e-5;
c_glucose = 5.55e-5*c;
c_water = 35.555; %35.555
c_lipid = 5.2e-5;

%% the extinction coefficient (unit: L*mol-1*cm-1)
sita_HGb = 1024;
sita_glucose = 1.48e4;
sita_water = 1.2e-3;
sita_lipid = 4.509;

output_args = log(10)*(c_HGb*sita_HGb + c_glucose*sita_glucose ...
    + c_water*sita_water + c_lipid*sita_lipid);
%% transform the unit to mm-1
output_args = output_args/10;
end

