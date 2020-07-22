%% Create a data set for square impacts
%number of trials you wish to collect data from
trials = 20;

%good trials... data andy and I collected from he visualizers manually
goodTrials = [find(trialData(:,1) == 1)',find(trialData(:,2) == 1)'];

%loop through trials
for k = 1:(length(goodTrials))
    %run the read data function to retrieve pre and post impact
    %positions/velocities as well as the normal and tangential Jacobians
    %for the impact
    
    i = goodTrials(k);
     
    [Pre, Post, d, n] = ReadData(i);
    squareDataPhil(k).states = [Pre, Post];
    squareDataPhil(k).n = n;
    squareDataPhil(k).d = d;
    squareDataPhil(k).trial = goodTrials(k);
end
