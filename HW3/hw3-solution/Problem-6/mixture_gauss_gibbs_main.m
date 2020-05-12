close all;
clear all;
clc; 

load mixture_gauss_basic_sample.mat; 


%% Parameters
burn_in_time = 10000; 
sample_number = 1000; 
sample_interval = 120; 

x_in = samples_record; 

samples = mix_gauss_unit_sample(x_in, burn_in_time, sample_number, sample_interval);