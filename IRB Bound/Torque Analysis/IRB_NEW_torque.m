
%{
This function computes the IRB Bound for a specified trial number
in order to compute the maximum area/width that the contact point could 
achieve. It takes into account the torque that is generated from the 
impact.
%}

function [P,w] = IRB_NEW_torque(trial)

load('ellipse_uniform.mat');

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
Mass = [m1, 0, 0;
        0, m1, 0; 
        0, 0, I1]; 


% Finding pre and post impact velocities / states
pre = bounce_array(trial).states(4:6)';
post = bounce_array(trial).states(10:12)';

d = (bounce_array(trial).d);   %tangential
n = (bounce_array(trial).n);   %normal

J = [d;n;0 0 1]; %Jacobian

fun = @(P)(findError_torque(P, Mass, J, pre, post));
nonlcon = @(P)(constraint_torque(P, Mass, J, pre, post));

P0 = [0 0 0];
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
P = fmincon(fun, P0, A, b, Aeq, beq, lb, ub, nonlcon);

w = abs(P(3))/P(2);

% disp("Tangential Impulse: " + P(1) + " [N*s]")
% disp("Normal Impulse: " + P(2) + " [N*s]")
% disp("Torque: " + P(3) + "[N*m]")
% disp("w: " + w + "[m]")

end