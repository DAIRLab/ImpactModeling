%% Find Best Mu and Epsilon using AP Poisson
clear;
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
sumX = [0, 0];
bests = [];
count = 0;
totalError = 0;
trialError = 0;
%set up mu and epsilon intervals for fmincon
mu_lim = [0, 0.4];
ep_lim = [0.2, 1];
lb = [mu_lim(1),ep_lim(1)];
ub = [mu_lim(2),ep_lim(2)];
%% Run Simulation
%loop over multiple trials
numTrials = 2000; %Number of Trials
ran = randi([1 2000], 1, numTrials);

for t = 3:numTrials
    trial = t;
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
    %Observed contact point momentum
    Pt = (M*(J*(post-pre)'));
    
    %if there are no flags run fmincon to minimize 
    if sum(bounce_array(trial).flags)<1
        %use Nima's fmincon settings 
        options = optimoptions('fmincon','FiniteDifferenceType','central',...
            'StepTolerance',1e-10,'Display','off');
        
        x0 = [0.1, 0.6]; %IC Nima uses after convexity check
        
        [x, fval] = fmincon(@(x)cost_function(x, Mass, M, n ,d, pre, post, J, Pt), x0, ...
            [], [], [], [], lb, ub, [], options);
        %keep track of # of entries
        count = count + 1;
        %add mu and epsilon to sum
        %sumX = sumX + x;
        trialError = cost_function(x, Mass, M, n ,d, pre, post, J, Pt);
        bests(count, :) = [x bounce_array(trial).states(3)];
        
        mu(count) = x(1);
        epsilon(count) =x(2);
        
        predicted = APPoisson_juniors(Mass, n, d, pre, min(mu), min(epsilon));
        
        errorVec(count) = norm(post - predicted');%/norm(post);
        
        useful(1, count) =  wrapTo180(rad2deg(bounce_array(trial).states(3)));
        useful(2, count) = norm(post);
    end
end

% figure()
% plot(useful(1,:), useful(2,:), '.');
% xlabel("Pre Impact Angle")
% ylabel("Post Impact Norm of Velocity");
% xlim([-180, 180]);
disp(mean(mu))
disp(mean(epsilon))
figure()
plot(useful(1,:), errorVec, '.');
xlabel("Pre Impact Angle [degrees]")
ylabel("l2 Norm Velocity Error");
xlim([-180, 180]);
%%
% 
% % find the avarage optimal mu and epsilon
% best = sumX/count;
% disp("Optimal Mu:  " + best(1))
% disp("Optimal Epsilon: " + best(2))
% 
% % %% Post Process Results
% % %Take the average of totalError
% averageError = totalError / count;
% minError = min(min(averageError));  %get the minimum error
% [a, b] = find(averageError == minError); %find its index
% 
%  %convert the index to mu and epsilon values
%  bestEp = sample(a);
%  bestMu = sample(b);
% % 
%  %create contour plot
%  colormap(flipud(gray))  %match color from Nima paper
%  contourf(sample(21:36), sample(1:15), averageError(21:36, 1:15)', 60);
%  colorbar; 
% % 
% xlabel("Epsilon")
% ylabel("Mu")
% 
% optTwos = 0;
% q = 0;
% optZeros = 0;
% r = 0;
% for h = 1:length(bests)
%     if bests(h, 1) > 0.05
%         optTwos = optTwos + abs(bests(h,3));
%         q = q + 1;
%     else
%         optZeros = optZeros + abs(bests(h,3));
%         r = r + 1;
%     end
% end
% 
% disp("Average observed w_dot for mu ~ 0.2: " + optTwos/q) 
% disp("Average observed w_dot for mu ~ 0: " + optZeros/r) 
