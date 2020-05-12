close all;
clear all;
clc; 

load mixture_gauss_basic_sample.mat; 

x_in = samples_record;
%% Parameters
burn_in_time = 10000; 
sample_number = 1000;
sample_interval = 30;

sigma_arr = [0.5 5] ;

for indx=1:2
    
    [samples, accept_rate] = mix_gauss_unit_metro_hast(x_in, burn_in_time, sample_number, sample_interval, sigma_arr(indx));
    fprintf('Sigma=%g,Acceptance rate=%g\n',sigma_arr(indx),accept_rate);

end

%% Result Documents 
% sigma = 0.5 mix_gauss_mp_s1_interval40
% Interval 10       20       30        40       50        60
% Accept   0.1271   0.1318   0.1291    0.1326   0.1310    0.1306


% sigma = 5.0  mix_gauss_mp_s2_interval40
% Interval 10       20       30        40       50        60
% Accept   0.0016   0.0013   0.0018    0.0015   0.0020    0.0017         
% For this scenario, it is more drastically changed. In some cases, we even
% have them to flip. 





























