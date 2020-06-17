%Run_explore_error

% This code takes n random trials (of 2000) and finds the 
% best mu and epsilon for each trial. After storing these
% in the vectors bmv and bev respectively (of size n), we 
% compute the mean and the standard deviation.

% Our first goal is to globally modify the pre
% positions, pre velocities, individually or at the same
% time, such that we get new bmv and bev vectors. To these, 
% we also compute the mean and standard deviation. 

% Our second goal is to try, as we change these variables,
% to obtain a standard deviation (for both bev and bmv) that 
% is each time lower - this means, our set of mus and epsilons
% for those n random trials are getting closer and closer to
% mean. This would be pretty ideal, since between two surfaces
% there should only be a single mu - hence, a lower standard
% deviation means we are getting closer to that

% Our third goal is, once we identify the parameters, to see if 
% they will reduce the standard deviation for all 2000 cases

% Step 1: set-up p1 and p2 (the max percentage change)
pmax1 = 0.5;
pmax2 = 0.5;
itr = 20; %how many iterations we would like
s1 = linspace(1-pmax1, 1+pmax1, itr); 
s2 = linspace(1-pmax2, 1+pmax2, itr); 
standard = zeros(itr,itr);
C = cell(itr,itr);

% Step 2: run the double for-loop
for i = 1:itr
    
    p1 = s1(i);

    %loop over all epsilons
    for j = 1:itr
        %assign mu to value corresponding with loop
        p2 = s2(j);
         
         % Calling on Exploring_Errors, which returns the average epsilon 
         % and mu for 1000 trials, as well as their standard deviations.
         % This is all data given the percentage change p1,p2 in the 
         % positions of impact
         [avEp,avMu,Standard_e,Standard_mu] = Exploring_Errors(p1,p2);
         standard(i, j) = Standard_mu;
         C{i,j} = {avEp,avMu,Standard_e,Standard_mu};
         
    end
    disp(i)
end

% Step 3: processing data
minStandard = min(min(standard));  %get the minimum standard
[a , b] = find(standard == minStandard); %find its index
disp(C{a,b})
p1 = s1(a)
p2 = s2(b)
