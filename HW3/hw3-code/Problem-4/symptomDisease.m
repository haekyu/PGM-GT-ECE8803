import brml.*

% Read SymptomDisease.mat and load b, p, s, W
load SymptomDisease.mat

% Get S: the number of sympotms and D: the number of diseases
[S, D] = size(W);

% Initialize Gibbs sampling
num_samples = 10000;
d = rand(D, 1) < p;
samples = zeros(D, num_samples);

% Gibbs sampling
for sample_i = 1 : num_samples
    % random reordering
    r = randperm(D); 

    for i = r
        d(i) = 0;
        tmp = sigmoid(W*d + b, 1);
        E0 = sum(d .* log(p) + (1-d) .* log(1-p)) + sum(s .* log(tmp)) + sum((1-s) .* log(1-tmp));

        d(i) = 1;
        tmp = sigmoid(W*d + b, 1);
        E1 = sum(d .* log(p) + (1-d) .* log(1-p)) + sum(s .* log(tmp)) + sum((1-s) .* log(1-tmp));
         
            % sample from the conditional p(d_i=1|,d_\i,s)
            if rand < sigmoid(E1-E0, 1) 
                d(i)=1;
            else 
                d(i)=0;
            end 
        
    end
    samples(:, sample_i) = d;
end

burnin=500;
mean(samples(:,burnin : end), 2)