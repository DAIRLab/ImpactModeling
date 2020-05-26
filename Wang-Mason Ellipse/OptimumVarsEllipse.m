% This code calulates the average best mu and epsilon for a number of
% trials. We can have the compiled code error function output if the
% trial is a sticking trial and the min mu for that specific trial. We 
% then have this code select the max mu for all sticking trials and weight
% that mu with the number of sticking trials, then go on and do an average
% of the mus and epsilons normally

function [OptMu, OptEps] = OptimumVarsEllipse
ran = randi([1 500],1,80);
MuVec = zeros(1,80);
EpsVec = zeros(1,80);
MuStickVec = [];
stickcount = 0;
for z = 1:80
    n = ran(z);
    [stick,Mu,Ep] = Error(n); %we can eedit output of this code
    if stick == 1
        MuStickVec(end+1) = Mu;
        sickcount = stickcount+1;
        EpsVec(z) = Ep;
    else
        MuVec(z) = Mu;
        EpsVec(z) = Ep;
    end
end
maxstickMu = max(MuStickVec);
MuStick_w = maxstickMu*stickcount; %weigh the value of the max mu of the sticking trials
[i,j] = size(MuVec);
OptMu = (sum(MuVec)+MuStick_w)/(j+stickcount);

[i,j] = size(EpsVec);
OptEps = sum(EpsVec)/j;
end