% This is the pseudo-code for the AP Poisson

load('ellipse_uniform.mat');

mu = 0.3;
e = 0.56;

tr = 1;
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
M = [m1,0,0;0,m1,0;0,0,I1];

pre = bounce_array(tr).states(4:6)';
d = (bounce_array(tr).d);
n = (bounce_array(tr).n);
%vf = vi + at - time here is 1/250
vf = pre + (1/250)*[0;-9.8;0];
eps = [1;1];
D = [-d;d];
Bnn = n*(M^-1)*n';
Bnt = n*(M^-1)*D';
Btn = Bnt';
Btt = D*(M^-1)*D';
bn = n*vf;
bt = D*vf;

A = [Bnn Bnt 0; Btn Btt eps; mu -eps' 0];
B = [Bnn Bnt; Btn Btt];
b = [bn; bt; 0]

[~,x] = LCPSolve(A,b);

%z = [Cn, beta1, beta2, lambda];

v = B*[x(1);x(2);x(3)] 
post = bounce_array(tr).states(10:12)'


