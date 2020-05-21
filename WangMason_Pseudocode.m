%Hi This is a starting pseudocode file to write down the flow of the code
%before we make changes to the main wang mason code file

%PS by initial quantities i mean preimpact quanitities


%define the constants B1, B2, B3: ALL are COM initial position (x1 and y1), masses, radius of gyration dependant.
% ^ equations: 19-21
%define the constants So and Co: both are initial COM velocities dependant (ydot_1o, ....) and initial COM position dependant (x1 and y1),
% ^ equations: 22-23


%figure out a way to find mu (e dependant)


% define Mu_s: B3 and B1 dependant
% ^ equation: 34

%define and identify small s (constant name s) from the inequalities using an if statement.
% ^ equation: 49
% define P_d and P_q:  Mu, the Bs, So and Co, and s dependant, all that are defined above
% ^equations: 35-36

%log in the measured post impact data
% define the contact point velocity (xdot_1c and ydot_1c): initial COM velocities dependant (ydot_1o, ....)
% ^paper 2 equations 8 and 9

% define the error metric in the first paper: dependant on the previous contact velocity 
% ^ first paper equations 4 to 7

%define e as a SYMBOL

%using poisson's model, find the contact mode: P_d, P_q, and e dependant: mu, mu_s, So, Co, s, e, and the Bs dependant
% ^ page 639 table 1

% ^ since we do not know e, we won't actually know what contact mode we are experiencing
% so, we define all the the resulting impulses Px and Py (in terms of e) for all contact modes and find for each contact mode
% the e that minimizes the metric error. and the contact mode with the e that results in the smallest metric error
%is the correct e with the most realistic contact mode.
% Since the Poisson Hypothesis of restitution says that the coefficient of restitution e is ONLY MATERIAL DEPENDANT, If we run
%this section multiple times for multiple trials (like what was done in paper 1, they "fit" the models to find the most
% realistic e and mu that fit 80 trials) we can find an e that is reasonable (or maybe find e values for all trials and then 
% average them?) 
% we can end this code here (aka a code to specifically find a good value of e) and then have another function that calls out
%this function to calculate the impulses and post impact state like mentioned below)


%so here after the previuos section, we have an e that minimizes the metric error on average for a good number of trials
%and then we can use that e to get the most realistic Px and Py given an input of preimpact velocities and positions.
% but all of this is mu dependant and I haven't thought of how to approach that yet.

%find the post impact velocities (xdot_1 and ...)
% ^ equations 5-7
