%IRB_bound

%Energy Constraint:
% 0.5*pre'*M*pre >= 0.5*(pre + (M^-1)*J'*P)'*M*(pre + (M^-1)*J'*P)

% Sizes:
% pre = 3x1 vector 
% M = 3x3 generalized mass matrix
% post = 3x1 vector
% J = 3x2 jacobian matrix
% P = 2x1 vector

%Goal: 
% 1) Find P such that we can create an energy ellipse diagram, 
% where the energy of the system is not exceeded
% 2) Minimize the error: min |post - (pre + (M^-1)*J'*P)|

% load in ellipse data 
load('ellipse_uniform.mat');

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
M = [m1,0,0;0,m1,0;0,0,I1]; 

% Setting up the trial number
trial = 5;

% Finding pre and post impact velocities
pre = bounce_array(trial).states(4:6)';
post = bounce_array(trial).states(10:12)';

% Get normal and tangetnial Jacobians from dataset
d = (bounce_array(trial).d);   %tangential
n = (bounce_array(trial).n);   %normal

J = [d;n];

% Establish initial energy
Initial = 0.5*pre'*M*pre;

% Applying fmincon
fun = @(P)norm(post - (pre + (M^-1)*J'*[P(1);P(2)]));
nonlcon = @(P) IRB(P, trial);
x0 = [0 0];
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
% output: x(1) - tangential impulse
%         x(2) - normal impulse
P = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)

%Plotting the energy ellipse
syms x y
assume(x >= 0)
assume(y >= 0)
p = 0.5*(pre + (M^-1)*J'*[x;y])'*M*(pre + (M^-1)*J'*[x;y]) - 0.5*pre'*M*pre;
a = expand(p)
ezplot(a, [-0.1 0.1 -0.02 0.23])
hold on 
plot([-0.1 0.1],[0 0],'b')
plot([0 0],[-0.02 0.23],'b')
title({'Energy Ellipse plot for trial #' num2str(trial)},'FontSize',15)
annotation('textbox',[.515 .7 .1 .2],'String','Py','EdgeColor','none','FontSize',15)
annotation('textbox',[.85 0.03 .1 .2],'String','Px','EdgeColor','none','FontSize',15)
hold on
plot(P(1),P(2),'*','MarkerSize',20,'linewidth',2)

% output: c - constraint (energy cannot increase)
function [c,ceq] = IRB(P, trial)
load('ellipse_uniform.mat');
% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
M = [m1,0,0;0,m1,0;0,0,I1]; 


% Finding pre and post impact velocities
pre = bounce_array(trial).states(4:6)';
post = bounce_array(trial).states(10:12)';

% Get normal and tangetnial Jacobians from dataset
d = (bounce_array(trial).d);   %tangential
n = (bounce_array(trial).n);   %normal

J = [d;n];

c(1) = 0.5*(pre + (M^-1)*J'*[P(1);P(2)])'*M*(pre + (M^-1)*J'*[P(1);P(2)]) - 0.5*pre'*M*pre;
ceq = [];
end
