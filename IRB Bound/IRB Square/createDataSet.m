%% Create a data set for square impacts
load('wedTD.mat');
%number of trials you wish to collect data from
trialData = wedTD;
%good trials... data andy and I collected from he visualizers manually
goodTrials = [find(trialData(:,1) == 1)', find(trialData(:,2) == 1)'];

%loop through trials
for k = 1:(length(find(trialData(:,1) == 1)))
    %run the read data function to retrieve pre and post impact
    %positions/velocities as well as the normal and tangential Jacobians
    %for the impact
    
    i = goodTrials(k);
    [Pre, Post, d, n] = ReadData(i, 1);
    squareDataPhilUpdated(k).states = [Pre, Post];
    squareDataPhilUpdated(k).n = n;
    squareDataPhilUpdated(k).d = d;
    squareDataPhilUpdated(k).trial = goodTrials(k);  
    squareDataPhilUpdated(k).impact = 1;
end


for j = 1:(length(find(trialData(:,2) == 1)))
    %run the read data function to retrieve pre and post impact
    %positions/velocities as well as the normal and tangential Jacobians
    %for the impact
    i = goodTrials(k+j);
    [Pre, Post, d, n] = ReadData(i, 2);
    squareDataPhilUpdated(k+j).states = [Pre, Post];
    squareDataPhilUpdated(k+j).n = n;
    squareDataPhilUpdated(k+j).d = d;
    squareDataPhilUpdated(k+j).trial = i;
    squareDataPhilUpdated(k+j).impact = 2;
end
