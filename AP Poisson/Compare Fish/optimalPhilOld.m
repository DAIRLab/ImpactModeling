%% Find Best Mu and Epsilon using AP Poisson

% load in ellipse data 
load('ellipse_uniform.mat');

%% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
Mass = [m1, 0, 0;
     0, m1, 0;
     0, 0, I1];

%% Set up interval
iter = 100; %how many iterations of mu and epsilon we would like
sample = linspace(0.01, 0.99, iter); %set up sample vector
trialError = zeros(iter);
totalError = zeros(iter);
count = 0;

bestMu = 0;
bestEp = 0;
%% Run Simulation
%loop over multiple trials
numTrials = 100; %Number of Trials
ran = randi([1 2000],1,numTrials);

for t = 1:numTrials
    trial = ran(t);
    %Get pre and post impact data from current trial
    pre = bounce_array(trial).states(4:6);
    post = bounce_array(trial).states(10:12);
    
    %get normal and tangetnial Jacobians from dataset
    d = (bounce_array(trial).d);   %tangential
    n = (bounce_array(trial).n);   %normal
    %variables for error calc
    J = [d;n];
    Minv = J*(Mass\J');
    M = inv(Minv);
    
    if sum(bounce_array(trial).flags)<1

        %loop over all mu's
        for i = 1:iter
            %assign epsilon to value corresponding with loop
            epsilon = sample(i);

            %loop over all epsilons
            for j = 1:iter
                %assign mu to value corresponding with loop
                mu = sample(j);

                %Run AP Poisson Model given mu and epsilon
                v_calc = APPoisson_juniors(Mass, n, d, pre, mu, epsilon);

                %calculate error for current trial using Nima's Method
                Pt = (M*(J*(post-pre)'));
                p_hat = (M*(J*(v_calc-pre')));
                error = norm(Pt - p_hat);
                trialError(i, j) = error;

            end
        end
        %add trial error to total error
        totalError = totalError + trialError;
        count = count + 1;
        
        minError = min(min(trialError));  %get the minimum error
        [a, b] = find(trialError == minError);
        bestMu = bestMu + sample(min(b));
        bestEp = bestEp + sample(min(a));
    end

end

%% Post Process Results
%Take the average of totalError
averageError = totalError / count;
minError = min(min(averageError));  %get the minimum error
[a, b] = find(averageError == minError); %find its index

disp(bestMu / count);
disp(bestEp / count);
%convert the index to mu and epsilon values
bestEp = sample(a);
bestMu = sample(b);

%create contour plot
colormap(flipud(gray))  %match color from Nima paper
contourf(sample(21:36), sample(1:15), averageError(21:36, 1:15)', 60);
colorbar; 

xlabel("Epsilon")
ylabel("Mu")