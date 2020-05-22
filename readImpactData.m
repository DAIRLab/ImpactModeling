%% Helper Functions

% This function takes in the specified .csv file and outputs the pre and
% post impact states of the first impact
% INPUTS: csvFile - name of the .csv file to be accessed
% OUTPUTS: pre - vector of pre impact state of first collision
%          post - vector of post impact state of first collision
function [pre, post] = readImpactData(csvFile) 

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