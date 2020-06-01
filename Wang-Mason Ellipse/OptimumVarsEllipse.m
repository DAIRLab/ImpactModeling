% This code calulates the average best mu and epsilon for a number of
% trials. We can have the compiled code error function output if the
% trial is a sticking trial and the min mu for that specific trial. We 
% then have this code select the max mu for all sticking trials and weight
% that mu with the number of sticking trials, then go on and do an average
% of the mus and epsilons normally

function [OptMu, OptEps] = OptimumVarsEllipse
ran = randi([1 2000],1,500);
MuVec = [];
EpsVec = [];
totalError = zeros(99,99);

% access actual data of first trajectory
load('ellipse_uniform.mat');

for n = 1:500
    z = ran(n);
    % pre - vector of pre impact state [x1_0, y1_0, theta1_0, x1dot_0, y1dot_0, theta1dot_0]
    % post - vector of post impact state [x1_act, y1_act, theta1_act, x1dot_act, y1dot_act, thetadot_act]
    pre = bounce_array(z).states(1:6); 
    post = bounce_array(z).states(7:12);
    J = [bounce_array(z).n; bounce_array(z).d];
    [stick,Mu,Ep,errors] = ErrorEllipse(z,pre,post,J); %we can eedit output of this code
    
    MuVec(end+1) = Mu;
    EpsVec(end+1) = Ep; 

    totalError = totalError + errors;
end

[i,j] = size(MuVec);
OptMu = mean(MuVec)

[i,j] = size(EpsVec);
OptEps = mean(EpsVec)

avgErrors = totalError./500;

mus = 0.01:0.01:0.99;
epsilons = 0.01:0.01:0.99;


contourf(epsilons, mus, avgErrors,20);
colorbar
hold on

%the best mu and epsilon for 500 random trials
scatter(EpsVec, MuVec, 15, 'k', 'filled')
hold on

%the optimal mu and epsilon overall (should be around (0.65, 0.16)
scatter(OptEps, OptMu, 20, 'w', 'filled')
hold off

xlabel('Epsilon')
ylabel('Mu')
title('Errors for Mu vs Epsilon')

end
