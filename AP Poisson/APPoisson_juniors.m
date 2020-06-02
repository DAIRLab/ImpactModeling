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


function [v_post] = APPoisson_juniors(M, n, d, v_pre, mu, epsilon)
%--- Set up variables ----
%number of edges of the polygonal approximation
u = 2; % two dimensional problem 
e = ones(2,1); %vectors of ones
D = [-d, d];










%Compression Phase

%Decompression Phase



end