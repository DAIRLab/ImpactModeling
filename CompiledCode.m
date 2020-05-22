%% Find Contact Mode and Apply Equations
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
stepSize = 0.05;
sl = 0.06; %side length of square from data README
rho = sqrt(sl^2/6); %using I/m where I = m *s^4 / 12
m1 = 1; %cancels out as explained in variables section above
I1 = m1 * sl^4 / 12; % moment of inertia of square
% initialize a matrix to hold the error values
errors = zeros(19,19); %size based on using 0.05 intervals
    
%% Find Impact Data
% access actual data of first trajectory

% pre - vector of pre impact x and y velocities [x1dot_0, y1dot_0]
% post - vector of post impact x and y velocities [x1dot_act, y1dot_act]
[pre, post] = actualVelocities('traj_1.csv'); 

%assign pre impact velocities for object 1
x1 = pre(1,1);
y1 = pre(1,2);
theta1 = pre(1,3);
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
    
for u = 0.05:stepSize:0.95 %varying mu from [0, 1] in intervals of 0.05
    Pd = (B2 + s * u * B3) * s * S_0;  %(35)
    Pq = (u * B1 + s * B3)*(-C_0);     %(36)
    
    for e = 0.05:stepSize:0.95 %varying epsilon from [0, 1] in intervals of 0.05
    
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
        x1dot_calc = Px/m1 + pre(1,4); 
        y1dot_calc = Py/m1 + pre(1,5); 
    
        % calculate error via least squares method ??
        error = (x1dot_calc - post(1,4))^2 + (y1dot_calc - post(1,5))^2; 
        % input error into error matrix
        errors(floor(u/stepSize), floor(e/stepSize)) = error;
        disp(floor(u/stepSize)+ " " + floor(e/stepSize)); 
    end
end

% determine minimum error 
    minimum = min(min(errors));
% determine the indices of the minimum error
    [i,j] = find(errors == minimum);
% determine the "best" mu and epsilon value which yields the minimum error
    bestMu = (i-1) * stepSize;
    bestEpsilon = (j-1) * stepSize;
    
%% Helper Functions

% This function takes in the specified .csv file and outputs the pre and
% post impact states of the first impact
% INPUTS: csvFile - name of the .csv file to be accessed
% OUTPUTS: pre - vector of pre impact state
%          post - vector of post impact state
function [pre, post] = actualVelocities(csvFile) 

    %load data
    traj = readtable(csvFile);
    
    %convert from table to matrix
    trajectory = traj{:, :};
    
    %initial impact matrix 
    impacts = zeros(1, 6, 2);

    %iterator variable to keep track of impact #
    curr = 1; 
    
    %loop through all of data 
    for i = 2:length(trajectory)-1
        if (trajectory(i-1,5) < 0) && (trajectory(i,5) > 0)
            impacts(curr, :, 1) = trajectory(i-1,:);  %preimpact data (first array)
            impacts(curr, :, 2) = trajectory(i,:);    %post impact data (second array)
            curr = curr + 1; 
        end
    end
    
    pre = impacts(1,:,1);
    post = impacts(1,:,2);
    
end
