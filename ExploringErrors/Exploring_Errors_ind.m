% Exploring Errors (individual version)

% This function uses the Wang model and finds the best mu and epsilon for
% a specific trial, given p1 and p2.

function [bestEp2,bestMu2,minError] = Exploring_Errors_ind(p1,p2,bounce_array,trial) 

% Step 1: set-up the constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
M = [m1,0,0;0,m1,0;0,0,I1];

% Step 2: set-up the interval
iter = 50; %how many iterations of mu and epsilon we would like
sample = linspace(0.01, 0.99, iter); %set up sample vector
trialError = zeros(iter);
totalError = zeros(iter);

% Step 3: run the simulation
%Get pre and post impact data from current trial
pre = bounce_array(trial).states(4:6);
post = bounce_array(trial).states(10:12);

%get normal and tangenial Jacobians from dataset
d = (bounce_array(trial).d);   %tangential
n = (bounce_array(trial).n);   %normal

J = [d;n];

%loop over all mu's
for i = 1:iter
    %assign epsilon to value corresponding with loop
    epsilon = sample(i);

    %loop over all epsilons
    for j = 1:iter
        %assign mu to value corresponding with loop
        mu = sample(j);

        %Run AP Poisson Model given mu and epsilon
         %v_calc = APPoisson_juniors(M, n, d, pre, mu, epsilon);
         v_calc = Wang_juniors(pre,n,d,mu,epsilon,p1,p2);

         error = sqrt((post(1) - v_calc(1))^2 + (post(2) - v_calc(2))^2)/sqrt(post(1)^2 + post(2)^2); 

         trialError(i, j) = error;            

    end

end


% Step 4: post-process results
%Finding average mu and ep by taking the average of the best of each trial
minError = min(min(trialError));
[be , bm] = find(trialError == min(min(trialError)));
bestEp2 = sample(be);
bestMu2 = sample(bm);

