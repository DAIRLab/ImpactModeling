%% Create a data set for square impacts
%load('updatedTD.mat');
%number of trials you wish to collect data from
trialData = updatedTD;
%good trials... data andy and I collected from he visualizers manually
goodTrials = [find(trialData(:,1) == 1)', find(trialData(:,2) == 1)'];

%loop through trials
for k = 1:(length(find(trialData(:,1) == 1)))
    %run the read data function to retrieve pre and post impact
    %positions/velocities as well as the normal and tangential Jacobians
    %for the impact
    
    i = goodTrials(k);
     
    [Pre, Post, d, n] = ReadData(i, 1);
    squareDataPhil(k).states = [Pre, Post];
    squareDataPhil(k).n = n;
    squareDataPhil(k).d = d;
    squareDataPhil(k).trial = goodTrials(k);
end

for k = 1:(length(find(trialData(:,1) == 1)))
    %run the read data function to retrieve pre and post impact
    %positions/velocities as well as the normal and tangential Jacobians
    %for the impact
    
    i = goodTrials(k);
     
    [Pre, Post, d, n] = ReadData(i, 1);
    squareDataPhil(k).states = [Pre, Post];
    squareDataPhil(k).n = n;
    squareDataPhil(k).d = d;
    squareDataPhil(k).trial = i;
    squareDataPhil(k).impact = 1;
end

for j = 1:(length(find(trialData(:,2) == 1)))
    %run the read data function to retrieve pre and post impact
    %positions/velocities as well as the normal and tangential Jacobians
    %for the impact
    
    i = goodTrials(k+j);
     
    [Pre, Post, d, n] = ReadData(i, 2);
    squareDataPhil(k+j).states = [Pre, Post];
    squareDataPhil(k+j).n = n;
    squareDataPhil(k+j).d = d;
    squareDataPhil(k+j).trial = i;
    squareDataPhil(k+j).impact = 2;
end
