% This code calulates the average best mu and epsilon for a number of
% trials. We can have the compiled code error function output if the
% trial is a sticking trial and the min mu for that specific trial. We 
% then have this code select the max mu for all sticking trials and weight
% that mu with the number of sticking trials, then go on and do an average
% of the mus and epsilons normally

function [OptMu, OptEps] = OptimumVarsEllipseFriday
ran = randi([1 2000],1,500);
MuVec = [];
EpsVec = [];
MuStickVec = [];
stickcount = 0;

% access actual data of first trajectory
load('ellipse_uniform.mat');

ErrorM = zeros(99,99);
minerr = [];

for z = 1:80
    n = ran(z);
    % pre - vector of pre impact state [x1_0, y1_0, theta1_0, x1dot_0, y1dot_0, theta1dot_0]
    % post - vector of post impact state [x1_act, y1_act, theta1_act, x1dot_act, y1dot_act, thetadot_act]
    pre = bounce_array(n).states(1:6); 
    post = bounce_array(n).states(7:12);
    J = [bounce_array(n).n; bounce_array(n).d];
    [stick,errors,Mu,Ep] = ErrorEllipseFriday(n,pre,post,J); %we can eedit output of this code
    
    if stick == 1
        MuStickVec(end+1) = Mu;
        stickcount = stickcount+1;
        EpsVec = [Ep,EpsVec];
    else
        MuVec = [Mu,MuVec];
        EpsVec = [Ep,EpsVec];
    end
    
    minimum = min(min(errors));
    
    minerr = [minimum,minerr];
    
    ErrorM = errors + ErrorM;
    
end
maxstickMu = max(MuStickVec);
MuStick_w = maxstickMu*stickcount; %weigh the value of the max mu of the sticking trials
[i,j] = size(MuVec);
OptMu = (sum(MuVec)+MuStick_w)/(j+stickcount)
UnweightedOptMu = mean([MuStickVec,MuVec])

[i,j] = size(EpsVec);
OptEps = sum(EpsVec)/j

ErrorM = ErrorM/z;
x = 0.01:0.01:.99;
y = 0.01:0.01:.99;
contourf(x,y,ErrorM,20)
colorbar

% minerr = zeros(1,99*99);
% 
% for k = 1:99
%     a = 99*k;
%     b = 99*(k-1)+1;
%     minerr(b:a) = ErrorM(k,:);
% end

%Creating the distribution plot
pd = fitdist(minerr','Normal');
step = 2/(length(minerr)-1);
x = 0:step:2;
y = pdf(pd,x);
plot(x,y,'LineWidth',2)
standardEp = pd.sigma;
yl = ylim;
%title({'Density Function for 80 Random Cases', 'Mean Error = ' num2str(avgerr)})
xlabel('Normalized Error Metric')
ylabel('Probability Density Estimate')

end
