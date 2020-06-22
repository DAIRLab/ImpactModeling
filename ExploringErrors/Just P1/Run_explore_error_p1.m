%Run_explore_error (only p1 version)

% This function goes trial by trial and computes the average mu and 
% as we change the x1 and y1 position by a percentage (p1 and p2). In 
% contrast to Run_explore_error where we assume a systematic error, and
% we try and find a p1 and p2 that reduce on average for ALL trials the
% standard deviation of mu, here we go individual trials - that way we can
% cater to each case, and find a distribution of the p1 and p2 

% Another big change is that now we are optimizing for the combination of
% p1 and p2 that yields the lowest error compared to finding the minimum
% standard deviation for mu

% Step 1: set-up p1 
pmax1 = 0.9;
itr = 19; %how many iterations we would like
s1 = linspace(1-pmax1, 1+pmax1, itr); 

% Step 2: load in ellipse data 
load('ellipse_uniform.mat');

% Step 3: set-up the number of trials
numTrials = 200;
errorM = zeros(1,itr,numTrials);
C = cell(1,itr,numTrials);
min_error_vec = zeros(3,numTrials);

% Step 4: run the double for-loop (outer: trials, inner: percentage
% changes in position)

for trials = 1:numTrials

    for i = 1:itr

        p1 = s1(i);
        p2 = 1;

        %loop over all values
         [avEp,avMu,minError] = Exploring_Errors_ind(p1,p2,bounce_array,trials);
         errorM(1, i,trials) = minError;
         C{1,i,trials} = {min(avEp),min(avMu),minError};
         
         if i == 10
           min_error_vec(2,trials) = minError;  
         end
        
    end
    disp(trials)
end

% Step 5: processing data (of all cases, regardless of whether they make
% physical sense or not)
percentage_vec = zeros(1,numTrials);
pos_vec = zeros(4,numTrials);
C2 = cell(numTrials,1);
for k = 1:numTrials  
    % Minimum error vector
    minEr = min(min(errorM(:,:,k)));  %get the minimum standard
    min_error_vec(1,k) = minEr;
    min_error_vec(3,k) = (minEr/min_error_vec(2,k))*100;
    % Percentage change vector
    a = find(errorM(:,:,k) == minEr); %find its index
    p1 = s1(a);
    percentage_vec(1,k) = p1;
    C2{k,1} = C{1,a,k};
    % Position vector
    d = (bounce_array(k).d);   %tangential
    n = (bounce_array(k).n); 
    J = [n;d];
    x1 = J(1,3);
    pos_vec(1:2,k) = [x1;x1*p1];
    pos_vec(3,k) = abs(pos_vec(1,k)-pos_vec(2,k));
    pos_vec(4,k) = (pos_vec(1,k)-pos_vec(2,k));
    end

meanOptimizedError = mean(min_error_vec(1,:))
meanRegularError = mean(min_error_vec(2,:))
avgPerChange = mean(min_error_vec(3,:))
mean_p1 = mean(percentage_vec(1,:))
meanPositionAbsolute = mean(pos_vec(3,:))
meanPositionRegular = mean(pos_vec(4,:))

% Step 6: making a scatter plot
plot(1000*pos_vec(3,:),min_error_vec(1,:),'o')
title(['Variation in error vs change in x1 position (absolute) for the first ' num2str(numTrials) ' trials'])
xlabel('Variation in position (mm)')
ylabel('Normalized Velocity error')
annotation('textbox',[.54 0.7 .4 .2],'String',['Avg. Actual error = ' num2str(meanRegularError)],'EdgeColor','none','FontSize',12)
annotation('textbox',[.54 0.665 .4 .2],'String',['Avg. Optimized error = ' num2str(meanOptimizedError)],'EdgeColor','none','FontSize',12)
annotation('textbox',[.54 0.63 .4 .2],'String',['Avg. Position (absolute) = ' num2str(meanPositionAbsolute)],'EdgeColor','none','FontSize',12)
set(gca, 'FontSize', 12)

figure
plot(1000*pos_vec(4,:),min_error_vec(1,:),'o')
title(['Variation in error vs change in x1 position (regular) for the first ' num2str(numTrials) ' trials'])
xlabel('Variation in position (mm)')
ylabel('Normalized Velocity error')
annotation('textbox',[.54 0.7 .4 .2],'String',['Avg. Actual error = ' num2str(meanRegularError)],'EdgeColor','none','FontSize',12)
annotation('textbox',[.54 0.665 .4 .2],'String',['Avg. Optimized error = ' num2str(meanOptimizedError)],'EdgeColor','none','FontSize',12)
annotation('textbox',[.54 0.63 .4 .2],'String',['Avg. Position (regular) = ' num2str(meanPositionRegular)],'EdgeColor','none','FontSize',12)
set(gca, 'FontSize', 12)

figure
plot(100*(1-percentage_vec(1,:)),min_error_vec(1,:),'o')
title(['Variation in error vs change (percentage) in x1 position for the first ' num2str(numTrials) ' trials'])
xlabel('Variation in position (%)')
ylabel('Normalized Velocity error')
annotation('textbox',[.57 0.7 .4 .2],'String',['Avg. Actual error = ' num2str(meanRegularError)],'EdgeColor','none','FontSize',12)
annotation('textbox',[.57 0.665 .4 .2],'String',['Avg. Optimized error = ' num2str(meanOptimizedError)],'EdgeColor','none','FontSize',12)
set(gca, 'FontSize', 12)

% figure
% h = histogram(1000*pos_vec(3,:),30);
% hold on
% plot(conv(h.BinEdges, [0.5 0.5], 'valid'), h.BinCounts)
