%Run_explore_error (individual version)
% Run time: ~7  minutes for 2000 trials

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
             C{i,j,trials} = {min(avEp),min(avMu),minError};

             if p1 == 1 & p2 == 1
               min_error_vec(2,trials) = minError;  
             end
             
        end
        
    end
    disp(trials)
end

% Step 5: processing data (of all cases, regardless of whether they make
% physical sense or not)
percentage_vec = zeros(2,numTrials);
pos_vec = zeros(4,numTrials);
C2 = cell(numTrials,1);
for k = 1:numTrials  
    % Minimum error vector
    minEr = min(min(errorM(:,:,k)));  %get the minimum standard
    min_error_vec(1,k) = minEr;
    min_error_vec(3,k) = (minEr/min_error_vec(2,k))*100;
    % Percentage change vector
    [a , b] = find(errorM(:,:,k) == minEr); %find its index
    p1 = s1(a);
    p2 = s2(b);
    percentage_vec(1:2,k) = [p1;p2];
    C2{k,1} = C{a,b,k};
    % Position vector
    d = (bounce_array(k).d);   %tangential
    n = (bounce_array(k).n); 
    J = [n;d];
    y1 = J(2,3);
    x1 = J(1,3);
    pos_vec(1:2,k) = [x1;x1*p1];
    pos_vec(3:4,k) = [y1;y1*p2];
    rad = sqrt((x1*p1)^2 + (y1*p2)^2);
    % Sanity Check - how many of our cases make physical sense? In other
    % words, does the new optimized x1 and y1 positions lie on the edge of
    % the ellipse (the radius from x1,y1 to the center of the ellipse is
    % between 0.025 and 0.0035)
    if (x1*p1 < 0.035) & (y1*p2 < 0.035) & (rad < 0.035 & rad > 0.025)
    pos_vec(5,k) = 1;
    else
    pos_vec(5,k) = 0;
    end
    
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
% 3) pos_vec: contains the values of x1, x1*p1, y1 and y1*p2 in  rows 1, 2,
% 3 and 4 respectively
% min_error_vec and percentage_vec have the same number of columns, which
% correspond to the trial number

% Observations: when running this code, almost for every trial the new
% optimized error < regular error, which means that just a slight change in
% positions (determined by p1 and p2) can make our model way more accurate.
% Furthermore, the values that are returned for p1 and for p2 are actually
% quite reasonable. This makes us wonder whether the data that we got or
% the approach to collecting it was actually well executed

% Step 6: post-processing of the physically sensible cases
vec = find(pos_vec(5,:) == 1);
l = length(vec);
percentage_vec2 = zeros(2,l);
pos_vec2 = zeros(5,l);
min_error_vec2 = zeros(3,l);
meanOptimizedError2 = mean(min_error_vec(1,vec));
meanRegularError2 = mean(min_error_vec(2,vec));
avgPerChange2 = mean(min_error_vec(3,vec));
mean2_p1 = mean(percentage_vec(1,vec));
mean2_p2 = mean(percentage_vec(2,vec));
percentage_vec2(:,:) = percentage_vec(:,vec);
pos_vec2(:,:) = pos_vec(:,vec);
min_error_vec2(:,:) = min_error_vec(:,vec);


% Note: when we isolate these cases, the statistics become even better -
% the average optimized error gets even smaller for instance.

% Another important thing to mention is that the # of isolated could
% potentially be expressed as a function of the percentage change. In other
% words, the more we let x1 and y1 vary (using p1 and p2), the less the #
% of isolated cases, since the pair won't be anymore within our desired
% area of focus
% For instance, for the first 200 cases, at a variation of p1 = 0.5 to 1.5
% and p2 = 0.5 to 1.5, we get 70 isolated cases (which is, the number of
% cases that make physical sense)
% In contast, if we limit p1 = 0.3 to 1.3
% and p2 = 0.3 to 1.3, we get 95 isolated cases. Also, at this range, we
% get a value for avgPerChange that is larger than before - hence, the more
% we can modify x1 and y1, the lower the errors we can get

% Some more statistics: p1 = 0.5 to 1.5 and p2 = 0.5 to 1.5, and total
% cases: 500, we get # isolated cases = 168, and avgPerChange = 24.2365
% (and 24.9757)

% When p1 = 0.5 to 1.5 and p2 = 0.5 to 1.5, and total
% cases: 1000, we get # isolated cases = 323, and avgPerChange = 24.0034
% (and 23.3808)

% When p1 = 0.5 to 1.5 and p2 = 0.5 to 1.5, and total
% cases: 2000, we get # isolated cases = 644, and avgPerChange = 23.2232
% (and 23.0018)

% Step 7: some graphs
plot(1:l,min_error_vec2(2,:))
hold on
plot(1:l,min_error_vec2(1,:))
legend('Actual error (p1, p2 = 1)','Optimized error (using p1, p2)')
title({'Actual vs Optimized error', ['Isolated Cases = ' num2str(l) ' of ' num2str(numTrials)]})
xlabel('Isolated Cases')
ylabel('Normalized Velocity error')
set(gca, 'FontSize', 12)
annotation('textbox',[.13 0.7 .4 .2],'String',['Avg. Actual error = ' num2str(meanRegularError2)],'EdgeColor','none','FontSize',12)
annotation('textbox',[.13 0.665 .4 .2],'String',['Avg. Optimized error = ' num2str(meanOptimizedError2)],'EdgeColor','none','FontSize',12)

% Step 8: mu/epsilon statistics for all cases
Stats = [C2{:,1}];
Stats = cell2mat(Stats);
AvgMu = mean(Stats(2:3:end));
AvgEp = mean(Stats(1:3:end));
% For the isolated cases
Stats2 = [C2{vec,1}];
Stats2 = cell2mat(Stats2);
AvgMu2 = mean(Stats2(2:3:end));
AvgEp2 = mean(Stats2(1:3:end));
