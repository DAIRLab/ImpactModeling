%% Find Best Mu and Epsilon using AP Poisson

% load in ellipse data 
load('ellipse_uniform.mat');

%% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
M = [m1, 0, 0;
     0, m1, 0;
     0, 0, I1];

%% Set up interval
iter = 50; %how many iterations of mu and epsilon we would like
sample = linspace(0.01, 0.99, iter); %set up sample vector
trialError = zeros(iter);
totalError = zeros(iter);


%% Run Simulation
%loop over multiple trials
numTrials = 1000; %Number of Trials
for trial = 1:numTrials
    %Get pre and post impact data from current trial
    pre = bounce_array(trial).states(4:6);
    post = bounce_array(trial).states(10:12);
    
    %get normal and tangetnial Jacobians from dataset
    d = (bounce_array(trial).d);   %tangential
    n = (bounce_array(trial).n);   %normal
    
    
    %loop over all mu's
    for i = 1:iter
        %assign epsilon to value corresponding with loop
        epsilon = sample(i);

        %loop over all epsilons
        for j = 1:iter
            %assign mu to value corresponding with loop
            mu = sample(j);

            %Run AP Poisson Model given mu and epsilon
            v_calc = APPoisson_juniors(M, n, d, pre, mu, epsilon);
            
            %calculate error for current trial
            error = sqrt(...%(post(1) - v_calc(1))^2 + ...
                (post(2) - v_calc(2))^2);%/sqrt(post(1)^2 + post(2)^2); 
            trialError(i, j) = error;
            
        end
    end
    %add trial error to total error
    totalError = totalError + trialError;
    
end

%% Post Process Results
%Take the average of totalError
averageError = totalError / numTrials;
minError = min(min(averageError));  %get the minimum error
[a, b] = find(averageError == minError); %find its index

%convert the index to mu and epsilon values
bestEp = sample(a);
bestMu = sample(b);

%create contour plot
colormap(flipud(gray))  %match color from Nima paper
contourf(sample(21:36), sample(1:15), averageError(21:36, 1:15)', 60);
colorbar; 

xlabel("Epsilon")
ylabel("Mu")




