% This code calulates the average best mu and epsilon for a number of
% trials. We can have the compiled code error function output if the
% trial is a sticking trial and the min mu for that specific trial. We
% then have this code select the max mu for all sticking trials and weight
% that mu with the number of sticking trials, then go on and do an average
% of the mus and epsilons normally

% function [OptMu, OptEps] = OptimumVarsEllipseFriday
ran = randi([1 2000],1,500);
MuVec = [];
EpsVec = [];
MuStickVec = [];
stickcount = 0;

% access actual data of first trajectory
load('ellipse_uniform.mat');
errorssum = zeros(99,99,1);

for n = 1:1000
    %n = ran(z);
    % pre - vector of pre impact state [x1_0, y1_0, theta1_0, x1dot_0, y1dot_0, theta1dot_0]
    % post - vector of post impact state [x1_act, y1_act, theta1_act, x1dot_act, y1dot_act, thetadot_act]
    pre = bounce_array(n).states(1:6);
    post = bounce_array(n).states(7:12);
    J = [bounce_array(n).n; bounce_array(n).d];
    [stick,Mu,Ep,errors] = ErrorEllipseFriday(n,pre,post,J); %we can eedit output of this code
    
    %     if stick == 1
    %         MuStickVec(end+1) = Mu;
    %         stickcount = stickcount+1;
    %         EpsVec = [Ep,EpsVec];
    %     else
    %         MuVec = [Mu,MuVec];
    %         EpsVec = [Ep,EpsVec];
    %     end
    errorssum = errorssum+ errors;
    
end

%getting the mu and the epsilon that minimize the eror sum
minimum = min(min(errorssum));
[i,j] = find(errorssum == minimum);
bMu = min(i * .01)
bEpsilon = min(j * .01)

%contour plot of the sum of the error matricies
muvec = 0.01:0.01:(1-0.01);
epvec = 0.01:0.01:(1-0.01);
contourf(epvec,muvec,errorssum,40)



%%getting the average mu and epsilon
% maxstickMu = max(MuStickVec);
% MuStick_w = maxstickMu*stickcount; %weigh the value of the max mu of the sticking trials
% [i,j] = size(MuVec);
% OptMu = (sum(MuVec)+MuStick_w)/(j+stickcount)
% UnweightedOptMu = mean([MuStickVec,MuVec])
%
% [i,j] = size(EpsVec);
% OptEps = sum(EpsVec)/j
% end