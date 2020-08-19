%% USING JOAH's METHOD, different version of createDataSet.m

%% Create a data set for square impacts
load('TD.mat');
%number of trials you wish to collect data from
trialData = TD;
%good trials... data andy and I collected from he visualizers manually
goodTrials = [find(trialData(:,1) == 1)', find(trialData(:,2) == 1)'];
count = 0;
%loop through trials
for k = 1:(length(find(trialData(:,1) == 1)))
    %run the read data function to retrieve pre and post impact
    %positions/velocities as well as the normal and tangential Jacobians
    %for the impact
    
    i = goodTrials(k);
    [Pre, Post, d, n] = findND(i, 1);
    if (n(3) ~= 0)
        count = count + 1;
        squareDataJoah(count).states = [Pre, Post];
        squareDataJoah(count).n = n;
        squareDataJoah(count).d = d;
        squareDataJoah(count).trial = i;  
        squareDataJoah(count).impact = 1;
    end
end


for j = 1:(length(find(trialData(:,2) == 1)))
    %run the read data function to retrieve pre and post impact
    %positions/velocities as well as the normal and tangential Jacobians
    %for the impact
    i = goodTrials(k+j);
    [Pre, Post, d, n] = findND(i, 2);
    if (n(3) ~= 0)
        count = count + 1;
        squareDataJoah(count).states = [Pre, Post];
        squareDataJoah(count).n = n;
        squareDataJoah(count).d = d;
        squareDataJoah(count).trial = i;
        squareDataJoah(count).impact = 2;
    end
end
