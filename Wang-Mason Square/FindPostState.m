% This function calculates the post-impact state of the square for a
% specified trial number, based on the best Mu and Epsilon values
% calculated from 80 random trials (calculated in OptimumVars.m)

function [x1dot_calc,y1dot_calc] = FindPostState(n)
%%Find Contact Mode and Apply Equations
%Authors: Joah, Phil, Natalie, Andy

% VARIABLES:
% stepSize - stepSize to vary the mu/epsilon by
% (x1dot_calc, y1dot_calc) - calculated post-impact x and y velocities
%                            of object 1
% (x1dot_act, y1dot_act) - actual post-impact x and y velocities of
%                          object 1 (from data set)
% (x1dot_0, y1dot_0, theta1dot_0) - pre-impact x, y, and angular
%                                   velocities of  object 1 (from data
%                                   set)
% (x2dot_0, y2dot_0, theta2dot_0) - pre-impact x, y, and angular
%                                   velocities of  object 2 which are
%                                   all always 0 becasue the table does
%                                   not move
% rho - radius of gyration of object 1
% I1 - moment of inertia of object 1
% S_0 - intial sliding velocity (calculated from data)
% C_0 - initial compression velocity (calculated from data)
% m1 - mass of object one, but since it will cancel itself out, we
%      simply set it equal to 1
%B1, B2, and B3 - constants dependant on geometry and mass properties
%                 of the system
% sl - side length of object being dropped
% s - initial sign of sliding velocity S0

%% Set up Variables
stepSize = 0.01;
sl = 0.06; %side length of square from data README
rho = sqrt(sl^2/6); %using I/m where I = m *s^4 / 12
m1 = 1; %cancels out as explained in variables section above
I1 = m1 * sl^2 / 6; % moment of inertia of square
% initialize a matrix to hold the error values
sz = 1/stepSize - 1;
errors = zeros(sz,sz); %size based on using 0.05 intervals

%% Find Impact Data
% access actual data of first trajectory

% pre - vector of pre impact x and y velocities [x1dot_0, y1dot_0]
% post - vector of post impact x and y velocities [x1dot_act, y1dot_act]
string = ['traj_' num2str(n) '.csv'];
[pre, post] = readImpactData(string);

%assign pre impact velocities for object 1
theta1 = pre(1,3);

%find the sign of x1, ie whether it is on the left or right of the contact point.
signx1 = 0;
if theta1>0 %turning anticlockwise
    r = rem(theta1,pi)/pi;
    if r >0 && r<0.25 || r >0.5 && r<0.75 %angle between 0 to 90 degrees and 180-270
        signx1 = 1; %x1 is to the right of contact point
    elseif r >0.25 && r<0.5 || r>0.75 %angle between 90 and 180 and 270-360
        signx1 = -1; %x1 will be on the left
    else %when it is right on top
        signx1 = 1;
    end
elseif theta1 <0 %turning clockwise
    r = rem(abs(theta1),pi)/pi;
    if r >0 && r<0.25 || r >0.5 && r<0.75 %angle between 0 to 90 degrees and 180-270
        signx1 = -1; %x1 is to the left of contact point
    elseif r >0.25 && r<0.5 || r>0.75 %angle between 90 and 180 and 270-360
        signx1 = 1; %x1 will be on the right
    else %when it is right on top
        signx1 = 1;
    end
elseif theta1 == 0 %doesn't even rotate
    signx1 = 1;
else
    signx1 = 1;
end

x1 = signx1*sqrt(sl^2/2 - pre(1,2)^2);
y1 = pre(1,2);
x1dot_0 = pre(1,4);
y1dot_0 = pre(1,5);
theta1dot_0 = pre(1,6);

%assign pre impact velocities for object 2
x2 = 1; %these are arbitrary since the velocities are 0 anyway
y2 = 1;
theta2 = 1;
x2dot_0 = 0;
y2dot_0 = 0;
theta2dot_0 = 0;

%% Solve for Constants
%S_0 and C_0
S_0 = (x1dot_0 + theta1dot_0 * y1) - (x2dot_0 + theta2dot_0 * y2); %(22)
C_0 = (y1dot_0 + theta1dot_0 * x1) - (y2dot_0 + theta2dot_0 * x2); %(23)
%B1, B2, and B3
B1 = 1 + y1^2/(m1*rho^2); %(19)
B2 = 1 + x1^2/(m1*rho^2); %(20)
B3 = x1 * y1/m1/rho^2;    %(21)
%Static friction coefficiant
u_s = -B3/B1; %(34)
%Initial Sign of Sliding Velocity
if (S_0 ~= 0)
    s = sign(S_0);
else
    s = 1;
end

%find optimum values of e and mu:
[u,e] = OptimumVars;

Pd = (B2 + s * u * B3) * s * S_0;  %(35)
Pq = (u * B1 + s * B3)*(-C_0);     %(36)

%Use Table 1 to determine modes (conditionals)
%Apply equations 39 - 48 based on mode
%Sliding (Second Row of Table)

if (Pd > (1+e)*Pq)
    Py = - (1+e) * C_0 / (B2 + s * u * B3); %(40)
    Px = - s * u * Py;                      %(39)
    
    %R (Third Row of Table)
elseif (Pq < Pd) && (Pd < (1+e) * Pq)
    if u > abs(u_s) %R-Sticking
        Py = -(1+e)* C_0 /(B2 + s * u * B3);  %(44)
        Px = (B3*Py - S_0) / B1;              %(43)
        
    else %R-Reversed Sliding
        Py = -(1+e) * C_0/(B2 + s * u * B3);              %(48)
        Px = s * u * (Py - 2 * S_0 / (B3 + s * u * B1));  %(47)
    end
    
    %C (Fourth Row of Table)
elseif (Pd < Pq)
    if u > abs(u_s) %C-Sticking
        Py = -(1+e) * (B1 * C_0 + B3 * S_0)/(B1*B2 - B3^2); %(42)
        Px = (B3 * Py - S_0)/B1;                            %(41)
    else %C-Reversed Sliding
        Py = -(1+e)/(B2 - s * u * B3)* (C_0 + ...
            (2*s*B3*S_0)/(B3 + s * u * B1));               %(46)
        Px = s * u *(Py - 2 * S_0 / (B3 + s * u * B1));    %(45)
    end
else
    disp("Error, none of contact mode requirements met");
end

% calculate the post impact velocities according to the contact
% mode
x1dot_calc = Px/m1 + pre(1,4)
y1dot_calc = Py/m1 + pre(1,5)

end