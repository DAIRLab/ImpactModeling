%% Find Impact Data
% access actual data of first trajectory

% pre - vector of pre impact x and y velocities [x1dot_0, y1dot_0]
% post - vector of post impact x and y velocities [x1dot_act, y1dot_act]
[pre, post] = actualVelocities('traj_2.csv'); 

%assign pre impact velocities for object 1

y1 = pre(1,2);
theta1 = pre(1,3);

omega = 2*pi-(2*abs(theta1)+(90*pi)/180);
xdistance = tan(abs(omega)/2)*y1;

if theta1 <0 %to the right of contact point
    x1 = abs(xdistance)
else %to the left
    x1 = - abs(xdistance)
end
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