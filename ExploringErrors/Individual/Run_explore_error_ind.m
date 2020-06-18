%Run_explore_error (individual version)

% This function goes trial by trial and computes the average mu and 
% as we change the x1 and y1 position by a percentage (p1 and p2). In 
% contrast to Run_explore_error where we assume a systematic error, and
% we try and find a p1 and p2 that reduce on average for ALL trials the
% standard deviation of mu, here we go individual trials - that way we can
% cater to each case, and find a distribution of the p1 and p2 

% Another big change is that now we are optimizing for the combination of
% p1 and p2 that yields the lowest error compared to finding the minimum
% standard deviation for mu

% Step 1: set-up p1 and p2 (the max percentage change)
pmax1 = 0.5;
pmax2 = 0.5;
itr = 11; %how many iterations we would like
s1 = linspace(1-pmax1, 1+pmax1, itr); 
s2 = linspace(1-pmax2, 1+pmax2, itr);

% Step 2: load in ellipse data 
load('ellipse_uniform.mat');

% Step 3: set-up the number of trials
numTrials = 100;
errorM = zeros(itr,itr,numTrials);
C = cell(itr,itr,numTrials);

% Step 4: run the triple for-loop (outer: trials, double-inner: percentage
% changes in position)

for trials = 1:numTrials

    for i = 1:itr

        p1 = s1(i);

        %loop over all epsilons
        for j = 1:itr
            %assign mu to value corresponding with loop
            p2 = s2(j);

             [avEp,avMu,minError] = Exploring_Errors_ind(p1,p2,bounce_array,trials);
             errorM(i, j,trials) = minError;
             C{i,j,trials} = {avEp,avMu,minError};

        end
        
    end
    disp(trials)
end

% Step 5: processing data
min_error_vec = zeros(1,numTrials);
percentage_vec = zeros(2,numTrials);
C2 = cell(numTrials,1);
for k = 1:numTrials  
minEr = min(min(errorM(:,:,k)));  %get the minimum standard
min_error_vec(k) = minEr;
[a , b] = find(errorM(:,:,k) == minEr); %find its index
p1 = s1(a);
p2 = s2(b);
percentage_vec(1:2,k) = [p1;p2];
c2{k,1} = C{a,b,k};
end

