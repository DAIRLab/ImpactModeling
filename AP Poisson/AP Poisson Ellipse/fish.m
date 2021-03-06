% This is the pseudo-code for the AP Poisson

load('ellipse_uniform.mat');

mu = 0.3; %coefficiant of static friction 
e = 0.56; %coefficiant of restitution   

tr = 190;

a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis

m1 = 0.037; %mass of ellipse [kg]

I1 = m1 * (a0^2 + b0^2) / 4; %moment of inertia (2D ellipse formula)
Mass = [m1,0,0;0,m1,0;0,0,I1]; %generalized mass matrix

ha = (1/250)*[0;-9.8;0]; %1/frequency * acceleration vector -- [m/s]

pre = bounce_array(tr).states(4:6)';

d = (bounce_array(tr).d)'; %tangential jacobian from data set
n = (bounce_array(tr).n)'; %normal jacobian from data set

u = 2; %number of edges of the polygonal approximation 

eps = [1;1]; %vector of ones same size as u

D = [-d,d]; %matrix with two columns, d and -d which represent vectors of 
            %tangent space
            
%Using Erleben's paper we can go from the Generalized LCP in Anitescu (eq. 
%13) to a true LCP where we can solve for the post impact state given the
%premimpact state and information about the geometry

%Set up LCP using equation 27 from Erleben (Matrix in LCP)
M = [n'*(Mass\n), n'*(Mass\D), 0;
     D'*(Mass\n), D'*(Mass\D), eps;
     mu, -eps', 0];

%(Vector in LCP)
%STILL NEED TO DOUBLE CHECK HOW WE GO FROM J * inv(M) * F to what we have
%below
q = [n'*(pre+ha); D'*(pre+ha); 0];

%Given matrix M and vector q, we now have a well defined LCP and can solve
%using an LCP solver
[w,z] = LCPSolve(M,q)