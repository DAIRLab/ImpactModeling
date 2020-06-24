%% Newer IRB

load('ellipse_uniform.mat');

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037;
I1 = m1 * (a0^2 + b0^2) / 4;
Mass = [m1, 0, 0;
    0, m1, 0;
    0, 0, I1];
Tlength = 101;
widthVec = zeros (1,Tlength);
TVec = zeros (1,Tlength);
PnVec = zeros (1,Tlength);
PtVec = zeros (1,Tlength);
errorMat = zeros (Tlength,101); %row trial, 101 width trials are columns
widthVecnonab = zeros (1,Tlength);


for i = 1:Tlength
    count = 0;
    for z = 0:0.1:10
        
        count= count+1
        %Select trial
        trial = i;
        
        width = z/1000;
        
        % Finding pre and post impact velocities / states
        pre = bounce_array(trial).states(4:6)';
        post = bounce_array(trial).states(10:12)';
        
        d = (bounce_array(trial).d);   %tangential
        n = (bounce_array(trial).n);   %normal
        
        J = [d;n;0 0 1]; %Jacobian
        
        fun = @(P)(findErrorwidth(P, Mass, J, pre, post, width));
        nonlcon = @(P)(constraintwidth(P, Mass, J, pre, post, width));
        
        P0 = [0 0]; %Pn Pt torque
        A = []; % No other constraints
        b = [];
        Aeq = [];
        beq = [];
        lb = [];
        ub = [];
        P = fmincon(fun, P0, A, b, Aeq, beq, lb, ub, nonlcon);
        
        errorMat(i,count) = findErrorwidth(P, Mass, J, pre, post, width);
        
        %     width = abs(P(3)) / P(2);
        %     widthnonab = P(3) / P(2);
        %
        %     widthVec(i) = width;
        %     widthVecnonab(i) = widthnonab;
        %
        %     errorVec(i) = findError(P, Mass, J, pre, post);
    end
end

% widthVec = widthVec*1000;
% avWid = mean(widthVec)
% avErr =  mean(errorVec)
AvErrVec = mean(errorMat);
widthVec = [0:0.1:10];
plot(widthVec,AvErrVec,".")


% disp("Tangential Impulse: " + P(1) + " [N*s]")
% disp("Normal Impluse: " + P(2) + " [N*s]")
% disp("Torque: " + P(3) + " [N*m]")
% disp("Average Width: " + avWid + " [mm]")
% disp("Average normalized error: " + avErr )
% histogram(errorVec,40)


