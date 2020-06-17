% Exploring Errors

function [bestEp2,bestMu2,Standard_e,Standard_mu] = Exploring_Errors(p1,p2,bounce_array) 

% Step 1: set-up the # of trials
numTrials = 1000; %Number of Trials

% Step 2: set-up the constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
M = [m1,0,0;0,m1,0;0,0,I1];

% Step 3: set-up the interval
iter = 50; %how many iterations of mu and epsilon we would like
sample = linspace(0.01, 0.99, iter); %set up sample vector
trialError = zeros(iter);
totalError = zeros(iter);

% Step 4: run the simulation
bev = zeros(1,numTrials);
bmv = zeros(1,numTrials);
%v = randi([1 2000],1,numTrials); %choose the random trials

for tr = 1:numTrials
    %trial = v(tr);
    trial = tr;
    %Get pre and post impact data from current trial
    pre = bounce_array(trial).states(4:6);
    post = bounce_array(trial).states(10:12);
    
    %get normal and tangetnial Jacobians from dataset
    d = (bounce_array(trial).d);   %tangential
    n = (bounce_array(trial).n);   %normal
    
    J = [d;n];
    Minv = J*(M\J');
    Ma = inv(Minv);
    
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
    %add trial error to total error
    totalError = totalError +  trialError;
    [be , bm] = find(trialError == min(min(trialError)));
    bev(tr) = min(be);
    bmv(tr) = min(bm);
end

% Step 5: post-process results
%Finding average mu and ep by taking the average of the best of each trial
bb = mean(bev);
aa = mean(bmv);
bestEp2 = sample(floor(bb)) + (sample(floor(bb)+1)-sample(floor(bb)))/(1)*(bb-floor(bb));
bestMu2 = sample(floor(aa)) + (sample(floor(aa)+1)-sample(floor(aa)))/(1)*(aa-floor(aa));

Standard_e = std(sample(bev));
Standard_mu = std(sample(bmv));
