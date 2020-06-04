% This is the pseudo-code for the AP Poisson

load('ellipse_uniform.mat');

mu = 0.3;
e = 0.56;

tr = 3;
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
M = [m1,0,0;0,m1,0;0,0,I1];

pre = bounce_array(tr).states(4:6)';
%Tangential Jacobian
d = (bounce_array(tr).d);
%Normal Jacobian
n = (bounce_array(tr).n);

%vf = vi + at - time here is 1/250
%Just like Prof. Posa said, think of this as an instantaneous moment where we 
%only care about velocities - the only acceleration we consider is gravity, 
%in the timespan of 1/250
v = pre + (1/250)*[0;-9.8;0];
eps = [1;1];
D = [-d;d];
Bnn = n*(M^-1)*n';
Bnt = n*(M^-1)*D';
Btn = Bnt';
Btt = D*(M^-1)*D';
%Note that here we multiply by n and D for the b vector to make sense (in other
%words, to have the correct size)
bn = n*v;
bt = D*v;

A = [Bnn Bnt 0; Btn Btt eps; mu -eps' 0];
B = [Bnn Bnt; Btn Btt];
b = [bn; bt; 0];

%Output x vector of normal and tangential impulses
[~,x] = LCPSolve(A,b);

%Poisson's Hypothesis - x(1) here corresponds to normal impulse
%x(2) and x(3) are the tangential impulses
x(1) = x(1) * (1 + e);

%Recall that mass*delta_v = Impulse

%This equation solves for delta_v and adds it to v - so we get the final v
%Note that we "divide" by n and D because earlier we had to multiply
v_calc = [M\n', M\D']*[x(1:3)]+(v)
post = bounce_array(tr).states(10:12)'
