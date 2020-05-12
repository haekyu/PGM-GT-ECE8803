close all;
clear all;
clc; 

sample_num = 100; 

samples_record = zeros(sample_num, 1);

mu = 5; 
sigma = 1; 

sample_func = @(mixture_idx) (1-mixture_idx).*normrnd(-1.0*mu, sigma, 1, 1) + ...
    mixture_idx .* normrnd(mu, sigma, 1, 1); 

for sample_idx = 1:1:sample_num
    % sample mixture idx 
    mix_idx = (rand(1) <= 0.5); 
    samples_record(sample_idx) = sample_func(mix_idx); 
end

hist(samples_record, 100)

save('mixture_gauss_basic_sample.mat', 'samples_record')

% !git add -A
% !git commit -m "mixture_gauss_basic_sample.mat"
% !git push