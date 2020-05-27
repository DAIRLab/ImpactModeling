ran = randi([1 2000],1,500);
MuVec = [];
EpsVec = zeros(1,500);
MuStickVec = [];
stickcount = 0;

% access actual data of first trajectory
load('ellipse_uniform.mat');


for z = 1:500
    n = ran(z);
    % pre - vector of pre impact state [x1_0, y1_0, theta1_0, x1dot_0, y1dot_0, theta1dot_0]
    % post - vector of post impact state [x1_act, y1_act, theta1_act, x1dot_act, y1dot_act, thetadot_act]
    pre = bounce_array(n).states(1:6); 
    post = bounce_array(n).states(7:12);
    [stick,Mu,Ep] = ErrorEllipse(n,pre,post); %we can eedit output of this code
    
    if stick == 1
        MuStickVec(end+1) = Mu;
        stickcount = stickcount+1;
        EpsVec(z) = Ep;
    else
        MuVec = [Mu,MuVec];
        EpsVec(z) = Ep;
    end
end
maxstickMu = max(MuStickVec);
MuStick_w = maxstickMu*stickcount; %weigh the value of the max mu of the sticking trials
[i,j] = size(MuVec);
OptMu = (sum(MuVec)+MuStick_w)/(j+stickcount)

[i,j] = size(EpsVec);
OptEps = sum(EpsVec)/j
end
