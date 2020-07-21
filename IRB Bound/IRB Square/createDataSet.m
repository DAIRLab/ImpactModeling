%% Create a data set for square impacts
%number of trials you wish to collect data from
trials = 20;

%loop through trials
for i = 1:trials
    %run the read data function to retrieve pre and post impact
    %positions/velocities as well as the normal and tangential Jacobians
    %for the impact
    [Pre, Post, d, n] = ReadData(i);
    squareData(i).states = [Pre, Post];
    squareData(i).n = n;
    squareData(i).d = d;
end
