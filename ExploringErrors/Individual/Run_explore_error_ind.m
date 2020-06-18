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
numTrials = 200;
errorM = zeros(itr,itr,numTrials);
C = cell(itr,itr,numTrials);
min_error_vec = zeros(3,numTrials);

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

             if p1 == 1 & p2 == 1
               min_error_vec(2,trials) = minError;  
             end
             
        end
        
    end
    disp(trials)
end

% Step 5: processing data
percentage_vec = zeros(2,numTrials);
C2 = cell(numTrials,1);
for k = 1:numTrials  
minEr = min(min(errorM(:,:,k)));  %get the minimum standard
min_error_vec(1,k) = minEr;
min_error_vec(3,k) = (minEr/min_error_vec(2,k))*100;
[a , b] = find(errorM(:,:,k) == minEr); %find its index
p1 = s1(a);
p2 = s2(b);
percentage_vec(1:2,k) = [p1;p2];
c2{k,1} = C{a,b,k};
end

meanOptimizedError = mean(min_error_vec(1,:))
meanRegularError = mean(min_error_vec(2,:))
avgPerChange = mean(min_error_vec(3,:))
mean_p1 = mean(percentage_vec(1,:))
mean_p2 = mean(percentage_vec(2,:))

% Output vectors:
% 1) min_error_vec: contains the optimized error after finding a p1 and p2
% that yields a value lower than the regular error (row 1), the regular
% error (row 2, when p1 and p2 both equal 1) and the percentage change,
% which is the optimized error/regular error (row 3)
% 2) percentage_vec: contains the values of p1 (row 1) and p2 (row 2) that
% yield the smallest error (as a reminder, p1 modifies x1 and p2 modifies
% x2)
% min_error_vec and percentage_vec have the same number of columns, which
% correspond to the trial number

% Observations: when running this code, almost for every trial the new
% optimized error < regular error, which means that just a slight change in
% positions (determined by p1 and p2) can make our model way more accurate.
% Furthermore, the values that are returned for p1 and for p2 are actually
% quite reasonable. This makes us wonder whether the data that we got or
% the approach to collecting it was actually well executed
