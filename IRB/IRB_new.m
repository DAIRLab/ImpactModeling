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

%Select trial
trial = 1;

% Finding pre and post impact velocities / states
pre = bounce_array(trial).states(4:6)';
post = bounce_array(trial).states(10:12)';

d = (bounce_array(trial).d);   %tangential
n = (bounce_array(trial).n);   %normal

J = [d;n]; %Jacobian

fun = @(x)(findError(x, Mass, J, pre, post));
nonlcon = @(x)(constraint(x, Mass, J, pre, post));

x0 = [0 0];
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
