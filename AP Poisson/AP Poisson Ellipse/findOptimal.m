%% Find Best Mu and Epsilon using AP Poisson

% load in ellipse data 
load('ellipse_uniform.mat');

%% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
M = [m1,0,0;0,m1,0;0,0,I1];
count = 0;
%% Set up interval
numTrials = 1;
ran = randi([1 2000],1,numTrials);
errorMatrix = zeros(99,99);
totalError = zeros(99,99);
xVel = zeros(99,99);
yVel = zeros(99,99);

for a = 1:numTrials
    z = ran(a);
    if sum(bounce_array(z).flags)<1

        % pre - vector of pre impact state [x1dot_0, y1dot_0, theta1dot_0]
        % post - vector of post impact state [x1dot_act, y1dot_act, thetadot_act]
        pre = bounce_array(z).states(4:6); 
        post = bounce_array(z).states(10:12);
        n = bounce_array(z).n; 
        d = bounce_array(z).d;

        for e = 1:99
            epsilon = e*0.01;

            for m = 1:99
                mu = m;

                %Run AP Poisson Model given mu and epsilon
                v_calc = APPoisson_juniors(M, n, d, pre, mu, epsilon);
                xVel(m, e) = v_calc(1);
                yVel(m, e) = v_calc(2);

                %calculate error
                error = findErrorFish(m, e, n, d, M, pre, post)
                errorMatrix(m, e) = error;
            end   
        end
        count = count + 1;
        [bmu, bep] = find(errorMatrix == min(errorMatrix));
        
        mu(count) = bmu;
        epsilon(count) = bep; 
        
    end
end

disp(mean(mu));
disp(mean(epsilon));

avgError = totalError / numTrials;
minError = min(min(avgError));

figure()
eps = 0.35:0.01:0.70;
mus = 0.01:0.01:0.30;
contourMatrix = avgError(1:30, 35:70);
contourf(eps, mus, contourMatrix, 30);
colorbar;
xlabel('Epsilon')
ylabel('Mu')
title('Error for Mu and Epsilon over 20 Random Trials')

[i, j] = find(avgError == minError);

optMu = i * 0.01
optEps = j * 0.01


