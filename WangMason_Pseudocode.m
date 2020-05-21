%%Overall Goal
%Step 1: Get initial velcoity values as well as information about object to solve for I
%Step 2: Solve for constants S0, C0, B1, B2, and B3 given values from step 1
%Step 3: Solve for post impact velocities given constants from step 2, mu, and epsilon
%Step 4: Perform a grid search over mu and epislon to minimize the error of the post impact velocities 
%Step 5: Return optimal mu, epsilon, and post impact velocities
hello
% ----------------------------------
%Hi This is a starting pseudocode file to write down the flow of the code
%before we make changes to the main wang mason code file

%PS by initial quantities i mean preimpact quanitities
%NOTE: Object 1 is a SQUARE with side length 0.06 meters

%Defining Variables:
    % (x1, y1) - COM position of object 1
    % p1 - Radius of Gyration of Object 1
    % I1 - Moment of Inertia of Object 1

%define the constants B1, B2, B3: ALL are COM initial position (x1 and y1), masses, radius of gyration dependant.
    B1 = (p1^2 + y1^2) / I1;
    B2 = (p1^2 + x1^2) / I1;
    B3 = x1 * y1 / I1;
% ^ equations: 19-21

%define the constants So and Co: both are initial COM velocities dependant (ydot_1o, ....) and initial COM position dependant (x1 and y1),
% ^ equations: 22-23

% define Mu_s: B3 and B1 dependant
    mu_s = -B3 / B1;
% ^ equation: 34

%define and identify small s (constant name s) from the inequalities using an if statement. only So dependant
% ^ equation: 49
% define P_d and P_q:  Mu, the Bs, So and Co, and s dependant, all that are defined above
% ^equations: 35-36

%For values of e from 0 to 1
    %for values of mu from 0 to 1
    
    %determine the contact mode using table 1 in the wang mason paper using if statements
    %determine the Px and Py for the determined contact mode
    
    %end
%end
    
%log in the measured post impact data
% define the contact point velocity (xdot_1c and ydot_1c): initial COM velocities dependant (ydot_1o, ....)
% ^paper 2 equations 8 and 9

% define the error metric in the first paper: dependant on the previously defined initial contact velocity 
% ^ first paper equations 4 to 7

%find the metric error value for each set of e and mu determined above for the givent trial 
%we do that by comparing the resulting values of Px and Py for each set with the measured Px and Py
%find the minimum error value within all sets and find the set of mu and e that gives that minimum error value
%repeat this for multiple trials to find optimal values of mu and e for a more general range of initial conditions
%then use this optimizing set (the set of mu and e) to calculate the value of Simulated Px and Py for needed trial

%find the post impact velocities (xdot_1 and ...) using the previously calculated Px and Py
% ^ equations 5-7
