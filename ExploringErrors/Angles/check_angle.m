
%load data
load('ellipse_uniform.mat');

%establish number of trials
numTrials = 2000;
check = zeros(1,numTrials);

for trial = 1:numTrials
    M = explore_angles(bounce_array,trial);
    if min(M(1,1:end-1)) < min(M(1,end))
        check(trial) = 1;
    else
        check(trial) = 0;
    end
    disp(trial)
end
        