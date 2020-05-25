% this code calulates the average best mu and epsilon for a number of
% trials. we can have the compiled code error function output the if the 
% trial is a sticking trial and the min mu for that specific trial. and
% then have this code select the max mu for all sticking trials and weight
% that mu with the number of stivking trials, then go on and do an average
% of the mus and epsilons normally 

ran = randi([1 500],1,80);
MuVec = zeros(1,80);
EpsVec = zeros(1,80);
for z = 1:80
    n = ran(z);
    [Mu,Ep] = Error(n); %we can eedit output of this code
    MuVec(z) = Mu;
    EpsVec(z) = Ep;
end
[i,j] = size(MuVec);
OptMu = sum(MuVec)/j;

[i,j] = size(EpsVec);
OptEps = sum(EpsVec)/j;