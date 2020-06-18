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

% Note: This code assumes there is a systematic mistake when
% reporting the values for position and velocity that is 
% repeated in every measurement

% Step 1: set-up p1 and p2 (the max percentage change)
pmax1 = 0.5;
pmax2 = 0.5;
itr = 20; %how many iterations we would like
s1 = linspace(1-pmax1, 1+pmax1, itr); 
s2 = linspace(1-pmax2, 1+pmax2, itr); 
standard = zeros(itr,itr);
C = cell(itr,itr);

% Step 2: load in ellipse data 
load('ellipse_uniform.mat');

% Step 3: run the double for-loop
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
         [avEp,avMu,Standard_e,Standard_mu] = Exploring_Errors(p1,p2,bounce_array);
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

% Results: turns out that for the first 1000 cases, using only Wang's model,
% we get that if the y position was 50% higher and the x position was 50% lower,
% then the standard deviation for the mu is the lowest.

% However, it should be noted that if we allowed the percentage changes to go
% higher than this (to 90% for instance), we would observe a similar pattern, in
% the sense that the standard deviation for mu is lowest when the y position is
% 90% higher and the x position is 90% lower.

% Furthermore, there is also an inverse relation between the standard deviation of
% epsilon and that of mu - the lower the standard for mu, the higher the one for e
% We are not really worrying too much about e because this can vary, since the 
% impact initial conditions dictates a lot about this
% In contrast, for mu, we should be expect a global value for it - this is why
% we are trying to reduce its standard deviation, since that would indicate
% most mus are close to each other.

% Statistics for first 1000 trials: code returns the cell {[0.4555] [0.0832] [0.1144]
% [0.0388]} where the first value represents average epsilon, then mu, and then the
% standard deviation for both e and mu - this value, of 0.0388 is the lowest that was
% reached for mu, again, assuming there was a systematic error with the way that the
% positions are measured
