function [samples, ret_accept_rate] = mix_gauss_unit_metro_hast(x_in, burn_in_time, sample_number, sample_interval, sigma)

mu_arr_old = zeros(2, 1);
mu_arr_new = zeros(2, 1);

samples = zeros(2, sample_number); 

accept_prob_const_func = @(mu1_old, mu2_old, mu1_new, mu2_new)...
    (mu1_old.^2 + mu2_old.^2 -  mu1_new.^2 - mu2_new.^2)./200;

accept_prob_func = @(val, mu1, mu2)...
    log(exp(-0.5.*(val-mu1).^2) + exp(-0.5.*(val-mu2).^2));

accept_prob = 0;

for bIdx = 1: 1: burn_in_time
    mu_arr_new(1) = normrnd(mu_arr_old(1), sigma);
    mu_arr_new(2) = normrnd(mu_arr_old(2), sigma); 
   
    %% Judge whether the jump is accepted or not.
    accept_prob = sum(accept_prob_func(x_in, mu_arr_new(1), mu_arr_new(2))) - ...
        sum(accept_prob_func(x_in, mu_arr_old(1), mu_arr_old(2)));
    accept_prob = exp(accept_prob + accept_prob_const_func(mu_arr_old(1), mu_arr_old(2), mu_arr_new(1), mu_arr_new(2)));
    
    accept_prob = min(1, accept_prob);
    
    if (rand(1) < accept_prob)
       mu_arr_old = mu_arr_new ;
    end
    
end

mu1_sum = 0; 
mu2_sum = 0; 
mu1_traj = zeros(sample_number, 1);
mu2_traj = zeros(sample_number, 1); 

ret_accept_rate = 0.0; 

for sIdx = 1:1:sample_number
   for iIdx = 1:1:sample_interval
    mu_arr_new(1) = normrnd(mu_arr_old(1), sigma);
    mu_arr_new(2) = normrnd(mu_arr_old(2), sigma); 
    
    %% Judge whether the jump is accepted or not.
    accept_prob = sum(accept_prob_func(x_in, mu_arr_new(1), mu_arr_new(2))) - ...
        sum(accept_prob_func(x_in, mu_arr_old(1), mu_arr_old(2)));
    accept_prob = exp(accept_prob + accept_prob_const_func(mu_arr_old(1), mu_arr_old(2), mu_arr_new(1), mu_arr_new(2)));
    
    accept_prob = min(1, accept_prob) ;
    
    if (rand(1) < accept_prob)
       mu_arr_old = mu_arr_new;  
       ret_accept_rate = ret_accept_rate + 1; 
    end
   end
   
   mu1_sum = mu1_sum + mu_arr_old(1); 
   mu2_sum = mu2_sum + mu_arr_old(2); 
   
   mu1_traj(sIdx) = mu1_sum/sIdx; 
   mu2_traj(sIdx) = mu2_sum/sIdx;
   
   samples(:, sIdx) = mu_arr_old;  
end

ret_accept_rate = ret_accept_rate/(sample_number * sample_interval);

figure; 
plot(mu1_traj, 'b'); 
hold on; 
plot(mu2_traj, 'r'); 
grid on;
xlabel('Sample Index');

end





