%% AP Poisson Model
%From Formulating Dynamic Multi-Rigid-Body Contact Problems with Friction 
% as Solvable Linear Complementarity Problems
% M. ANITESCU
% F. A. POTRA

%This problems solves the LCP to resolve post impact states given the
%preimpact states as well as information about the objects geometry

% Inputs
% M - mass matrix (3x3 Identity matrix but with [m, m, I] on the diagonals) 
% d - tangential Jacobian (from data set)
% n - normal Jacobian (from data set) 
% v_pre - intial velocity of body including xdot, ydot, and thetadot
% mu - coefficiant of friction (Coulumbic Friction)
% epsilon - coefficiant of resitution


%Outpts
%v_post - the post impact velocity 


function [v_post] = APPoisson_juniors(M, n, d, pre, mu, epsilon)
    %--- Set up variables ----
    %number of edges of the polygonal approximation
    vf = pre' + (1/250)*[0;-9.8;0];
    eps = [1;1];
    D = [-d;d];
    Bnn = n*(M^-1)*n';
    Bnt = n*(M^-1)*D';
    Btn = Bnt';
    Btt = D*(M^-1)*D';
    bn = n*vf;
    bt = D*vf;

    A = [Bnn Bnt 0; Btn Btt eps; mu -eps' 0];
    %B = [Bnn Bnt; Btn Btt];
    b = [bn; bt; 0];

    %Output x vector of normal and tangential impulses
    [~,x] = LCPSolve(A,b);

    %Poisson's Hypothesis - x(1) here corresponds to normal impulse
    %x(2) and x(3) are the tangential impulses
    x(1) = x(1) * (1 + epsilon);

    %Recall that mass*delta_v = Impulse

    %This equation solves for delta_v and adds it to v - so we get the final v
    %Note that we "divide" by n and D because earlier we had to multiply
    v_post = [M\n', M\D']*[x(1:3)]+(vf);
    

end