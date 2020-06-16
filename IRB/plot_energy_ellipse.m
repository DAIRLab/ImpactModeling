% Plots the energy ellipse with the Pt and Pn from fmincon
function plot_energy_ellipse(trial)

syms x y
assume(x >= 0)
assume(y >= 0)

load('ellipse_uniform.mat');

% %run IRB script to find optimal impulse
% run('IRB_bound.m');

impulse = IRB_NEW(trial);

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
M = [m1,0,0;0,m1,0;0,0,I1]; 

% Setting up the trial number
% trial = 20;

% Finding pre and post impact velocities
pre = bounce_array(trial).states(4:6)';
post = bounce_array(trial).states(10:12)';

% Get normal and tangetnial Jacobians from dataset
d = (bounce_array(trial).d);   %tangential
n = (bounce_array(trial).n);   %normal

J = [d;n];

p = 0.5*(pre + (M^-1)*J'*[x;y])'*M*(pre + (M^-1)*J'*[x;y]) - 0.5*pre'*M*pre;
a = expand(p)
ezplot(a, [-0.1 0.1 -0.02 0.23])
hold on 
plot([-0.1 0.1],[0 0],'b')
plot([0 0],[-0.02 0.23],'b')
hold on
scatter(impulse(1), impulse(2), 'k', 'filled')
title({'Energy Ellipse plot for trial #' num2str(trial)},'FontSize',15)
annotation('textbox',[.515 .7 .1 .2],'String','Py','EdgeColor','none','FontSize',15)
annotation('textbox',[.85 0.03 .1 .2],'String','Px','EdgeColor','none','FontSize',15)

end
