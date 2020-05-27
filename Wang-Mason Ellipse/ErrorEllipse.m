% This function calculates the best Mu and best Epsilon for one trial


function [stick, bestMu, bestEpsilon] = ErrorEllipse(n)

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
    % B1, B2, and B3 - constants dependant on geometry and mass properties
    %                 of the system
    % sl - side length of object being dropped 
    % s - initial sign of sliding velocity S0
    
%% Set up Variables
stepSize = 0.01;
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
rho = 0.5 * sqrt(a0^2 + b0^2); %using rho = sqrt(I/m) where I = m*(a^2 + b^2)/4
m1 = 1; %cancels out as explained in variables section above
I1 = m1 * (a0^2 + b0^2) / 4; % moment of inertia of elliptical disk
% initialize a matrix to hold the error values
sz = 1/stepSize - 1;
errors = zeros(sz,sz); %size based on using intervals that exclude 0 and 1
    
%% Find Impact Data
% access actual data of first trajectory
load('ellpse_uniform.mat');

% pre - vector of pre impact state [x1_0, y1_0, theta1_0, x1dot_0, y1dot_0, theta1dot_0]
% post - vector of post impact state [x1_act, y1_act, theta1_act, x1dot_act, y1dot_act, thetadot_act]
[pre] = bounce_array(n).states(1:6); 
[post] = bounce_array(n).states(7:12);

%assign pre impact velocities for object 1
theta1 = pre(1,3);
y1 = pre(1,2); 


%find the sign of x1, ie whether it is on the left or right of the contact point.
signx1 = 0;
if theta1 > 0 %turning anticlockwise
    r = rem(theta1,pi)/pi;
    if r > 0 && r < 0.25 || r > 0.5 && r < 0.75 %angle between 0 to 90 degrees and 180-270
        signx1 = 1; %x1 is to the right of contact point
    elseif r > 0.25 && r < 0.5 || r > 0.75 %angle between 90 and 180 and 270-360
        signx1 = -1; %x1 will be on the left
    else %when it is right on top
        signx1 = 1;
    end
elseif theta1 < 0 %turning clockwise
    r = rem(abs(theta1),pi)/pi;
    if r > 0 && r < 0.25 || r > 0.5 && r < 0.75 %angle between 0 to 90 degrees and 180-270
        signx1 = -1; %x1 is to the left of contact point
    elseif r > 0.25 && r < 0.5 || r > 0.75 %angle between 90 and 180 and 270-360
        signx1 = 1; %x1 will be on the right
    else %when it is right on top
        signx1 = 1;
    end
elseif theta1 == 0 %doesn't even rotate
    signx1 = 1;
else
    signx1 = 1; 
end

%use a tilted ellipse equation to find x1 (x distance between the COM and the contact point)

syms x y h k a b th

eq1 = (((x-h)*cos(th) + (y-k)*sin(th))^2)/a^2 + (((x-h)*sin(th) - (y-k)*cos(th))^2)/b^2 == 1;
x1eq = solve(eq1, sym('x'));


a0 = 0.07/2; %semi axis major
b0 = 0.05/2; %semi axis minor
h0 = 0; %h0 and k0 are 0 since we assume the com is the origin of the ellipse
k0 = 0;
th0 = theta1; %tilting angle
y0= -y1; %distance from com to contact point (going down from com to contact point hence the negative sign) because y1 
%goes up from the table ie contact point to COM (

%solving for x1
x1 = subs(x1eq, [a,b,h,k,th,y], [a0,b0,h0,k0,th0,y0]); %gives out the x distance from the COM to the contact point in both positive
%negative since the answer is a square root (a vector of 2 elements)

x1 = signx1*abs(x1(1)); %since we figured out the sign of x1 above in the if statement, we multiply it by the the absolute of the 
%value of the first elemtent of the x1 vector. (it doesn't matter whether we choose the first or second elements because they
%have the same value just different signs


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


%% use error metric to determine best mu and e
% 1. compute error of the matrix containing Px & Py values by comparing to actual data
% 2. determine the minimum error and its position in the matrix
% 3. based on the position in the matrix, determine the mu and e values

% ASSUMPTIONS:
    % we will use a matrix containing calculated Px and Py values
    % columns will be varying epsilon
    % rows will be varying mu
    
for a = 1:sz %varying mu from [0, 1]
    u = stepSize * a;
    
    Pd = (B2 + s * u * B3) * s * S_0;  %(35)
    Pq = (u * B1 + s * B3)*(-C_0);     %(36)
    
    for b = 1:sz %varying epsilon from [0, 1] in intervals of 0.05.05, 0.95]
        %Use Table 1 to determine modes (conditionals)
        %Apply equations 39 - 48 based on mode
        %Sliding (Second Row of Table)
        e = stepSize * b;
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
        x1dot_calc = Px/m1 + pre(1,4); 
        y1dot_calc = Py/m1 + pre(1,5); 
    
        % calculate error via least squares method ??
        error = ((x1dot_calc - post(1,4))^2 + (y1dot_calc - post(1,5))^2)/(post(1,4)^2 + post(1,5)^2); 
        % input error into error matrix
        errors(a, b) = error;
    end
end

% determine minimum error 
    minimum = min(min(errors));
% determine the indices of the minimum error
    [i,j] = find(errors == minimum);
% determine the "best" mu and epsilon value which yields the minimum error
    bMu = i * stepSize;
    bEpsilon = j * stepSize;
    
    [g,h] = size(bMu);
    if g>1 || h>1 
        stick = 1; %multiple best mus - sticking
    else
        stick = 0; %a single best my - no sticking
    end
    
    bestMu = min(bMu);
    bestEpsilon = min(bEpsilon);
    
end
