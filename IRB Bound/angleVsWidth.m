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
Tlength = 1000;
widthVec = zeros (1,Tlength);
% TVec = zeros (1,Tlength);
% PnVec = zeros (1,Tlength);
% PtVec = zeros (1,Tlength);
% errorVec = zeros (1,Tlength);
widthVecab = zeros (1,Tlength);
thetaVec = zeros(1,Tlength);

    
for i = 1:Tlength
    
%Select trial
trial = i;

% Finding pre and post impact velocities / states
pre = bounce_array(trial).states(4:6)';
post = bounce_array(trial).states(10:12)';

%theta
theta = bounce_array(trial).states(3);

if abs(theta)> pi
    theta = abs((rem(theta,pi)*180)/pi)
else
    theta = abs(theta*180/pi)
end

thetaVec(trial) = theta;


d = (bounce_array(trial).d);   %tangential
n = (bounce_array(trial).n);   %normal

J = [d;n; 0 0 1]; %Jacobian

fun = @(P)(findError(P, Mass, J, pre, post));
nonlcon = @(P)(constraint(P, Mass, J, pre, post));

P0 = [0 0 0]; %Pn Pt torque
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
P = fmincon(fun, P0, A, b, Aeq, beq, lb, ub, nonlcon);

width = P(3) / P(2);
widthVec(i) = width;

widthab = abs(P(3)) / P(2);
widthVecab(i) = widthab;

errorVec(i) = findError(P, Mass, J, pre, post);

end

widthVec = widthVec*1000;
avWid = mean(widthVec)
avErr =  mean(errorVec)


% disp("Tangential Impulse: " + P(1) + " [N*s]")
% disp("Normal Impluse: " + P(2) + " [N*s]")
% disp("Torque: " + P(3) + " [N*m]")
% disp("Average Width: " + avWid + " [mm]")
% disp("Average normalized error: " + avErr )
% histogram(errorVec,40)

figure 1
plot(thetaVec,widthVec,".");
ylabel("Width / mm");
xlabel("theta / degree");

figure 2
plot(thetaVec,abs(widthVec),".");
ylabel("absolute Width / mm");
xlabel("theta / degree");

