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
