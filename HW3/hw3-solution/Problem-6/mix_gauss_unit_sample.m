function [samples] = mix_gauss_unit_sample(x, burn_in_time, sample_number, sample_interval)

z = zeros(length(x), 1); 

% initialize mu 
mu_arr = normrnd(0, 1, [1, 2]); 

% initialize z
z_func = @(val, mu_1, mu_2) ...
    2.0 - (rand(1) < exp(-0.5.*(val - mu_1).^2)./(exp(-0.5.*(val - mu_1).^2) + exp(-0.5.*(val - mu_2).^2)));

z_arr = z_func(x, mu_arr(1), mu_arr(2));

samples = zeros(length(x) + 2, sample_number); 


%% Define Parameters
mu1_traj = zeros(sample_number, 1); 
mu2_traj = zeros(sample_number, 1); 

mu1_sum = 0; 
mu2_sum = 0; 


%% burn_in time 

for burnIdx = 1:1:burn_in_time
% Update mu 

A_1 = 0.01 + sum(z_arr == 1); 
A_2 = 0.01 + sum(z_arr == 2); 

B_1 = sum(x(z_arr == 1)); 
B_2 = sum(x(z_arr == 2));

mu_arr(1) = normrnd(B_1./A_1, 1.0/sqrt(A_1));
mu_arr(2) = normrnd(B_2./A_2, 1.0/sqrt(A_2));

% update z
z_arr = z_func(x, mu_arr(1), mu_arr(2)); 

end

%% Pick samples 
for sIdx = 1:1:sample_number
   for iIdx = 1:1:sample_interval
       
        A_1 = 0.01 + sum(z_arr == 1); 
        A_2 = 0.01 + sum(z_arr == 2); 

        B_1 = sum(x(z_arr == 1)); 
        B_2 = sum(x(z_arr == 2));

        mu_arr(1) = normrnd(B_1./A_1, 1.0/sqrt(A_1));
        mu_arr(2) = normrnd(B_2./A_2, 1.0/sqrt(A_2));

        % update z
        z_arr = z_func(x, mu_arr(1), mu_arr(2)); 
   end
   samples(:, sIdx) = [mu_arr.'; z_arr]; 
   
   mu1_sum = mu1_sum + mu_arr(1); 
   mu2_sum = mu2_sum + mu_arr(2); 
   
   mu1_traj(sIdx) =  mu1_sum/sIdx;  
   mu2_traj(sIdx) =  mu2_sum/sIdx; 
end

plot(mu1_traj, 'r');
hold on;
plot(mu2_traj, 'b');
xlabel('Sample Index')






end



















