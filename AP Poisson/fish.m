% This is the pseudo-code for the AP Poisson

load('ellipse_uniform.mat');

mu = 0.3;
e = 0.56;

tr = 190;
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
Mass = [m1,0,0;0,m1,0;0,0,I1];
ha = (1/250)*[0;-9.8;0];

pre = bounce_array(tr).states(4:6)';
d = (bounce_array(tr).d)';
n = (bounce_array(tr).n)';

eps = [1;1];
D = [-d,d];

M = [n'*(Mass\n), n'*(Mass\D), 0;
     D'*(Mass\n), D'*(Mass\D), eps;
     mu, -eps', 0];

q = [n'*(pre+ha); D'*(pre+ha); 0];

[w,z] = LCPSolve(M,q)