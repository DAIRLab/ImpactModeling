%Run_explore_error (Px/Py version)

% This function goes trial by trial and computes the average mu and 
% as we change the x1 position by a percentage (p1). In 
% contrast to Run_explore_error where we assume a systematic error, and
% we try and find a p1 that reduces on average for ALL trials the
% standard deviation of mu, here we go individual trials - that way we can
% cater to each case, and find a distribution of the p1  

% Step 1: set-up p1 
pmax1 = 0.9;
itr = 19; %how many iterations we would like
s1 = linspace(1-pmax1, 1+pmax1, itr); 

% Step 2: load in ellipse data 
load('ellipse_uniform.mat');

% Step 3: set-up the number of trials
desiredTrials = 1000;
vec = zeros(1,desiredTrials);
% double for loop to discard the flagged trials
for h = 1:desiredTrials
    if sum(bounce_array(h).flags) == 0
        vec(h) = 1;
    else
        vec(h) = 0;
    end
end
        
nonflagged = find(vec == 1);
numTrials = sum(vec == 1);
errorM = zeros(1,itr,numTrials);
errorP = zeros(2,itr,numTrials);
C = cell(1,itr,numTrials);
min_error_vec = zeros(3,numTrials);

% Step 4: run the double for-loop (outer: trials, inner: percentage
% changes in position)

for tr = 1:numTrials
    trials = nonflagged(tr);
    for i = 1:itr

        p1 = s1(i);
        p2 = 1;

        %loop over all values
         [avEp,avMu,minError,P] = Exploring_Errors_p(p1,p2,bounce_array,trials);
         errorM(1, i,tr) = minError;
         errorP(:, i,tr) = P';
         C{1,i,tr} = {min(avEp),min(avMu),minError};

         if i == 10
           min_error_vec(2,tr) = minError;  
         end
    end
        
    disp(tr)
end

% Step 5: processing data (of all cases, regardless of whether they make
% physical sense or not)
percentage_vec = zeros(1,numTrials);
pos_vec = zeros(8,numTrials);
impulses = zeros(2,numTrials);
C2 = cell(numTrials,1);

for k = 1:numTrials  
    % Minimum error vector
    minEr = min(min(errorM(:,:,k)));  %get the minimum standard
    min_error_vec(1,k) = minEr;
    min_error_vec(3,k) = (minEr/min_error_vec(2,k))*100;
    % Percentage change vector
    a = min(find(errorM(:,:,k) == minEr)); %find its index
    p1 = s1(a);
    percentage_vec(1,k) = p1;
    impulses(1:2,k) = errorP(1:2,a,k);
    C2{k,1} = C{1,a,k};
    % Position vector
    d = (bounce_array(k).d);   %tangential
    n = (bounce_array(k).n); 
    J = [n;d];
    x1 = J(1,3);
    pos_vec(1:2,k) = [x1;x1*p1];
    pos_vec(3,k) = abs(pos_vec(1,k)-pos_vec(2,k));
    pos_vec(4,k) = (pos_vec(1,k)-pos_vec(2,k));
    % Y positions
    y0 = bounce_array(nonflagged(k)).states(2);
    yc = d(3);
    pos_vec(5,k) = yc;
    post_y  = bounce_array(nonflagged(k)).states(8);
    pos_vec(6,k) = post_y;
    % post y position
    post_x = bounce_array(nonflagged(k)).states(7);
    pos_vec(7,k) = post_x-x1*p1;
    % pre angle impact
    angle = bounce_array(nonflagged(k)).states(3);
    pos_vec(8,k) = angle;
    end

meanOptimizedError = mean(min_error_vec(1,:))
meanRegularError = mean(min_error_vec(2,:))
avgPerChange = mean(min_error_vec(3,:))
mean_p1 = mean(percentage_vec(1,:))
meanPositionAbsolute = mean(pos_vec(3,:))*1000
meanPositionRegular = mean(pos_vec(4,:))*1000

% Step 6: making a scatter plot
%Error vs absolute change in x1
plot(1000*pos_vec(3,:),min_error_vec(1,:),'o')
title(['Variation in error vs change in x1 position (absolute) for the first ' num2str(numTrials) ' unflagged trials'])
xlabel('Variation in position (mm)')
ylabel('Normalized Velocity error')
annotation('textbox',[.54 0.7 .4 .2],'String',['Avg. Actual error = ' num2str(meanRegularError)],'EdgeColor','none','FontSize',12)
annotation('textbox',[.54 0.665 .4 .2],'String',['Avg. Optimized error = ' num2str(meanOptimizedError)],'EdgeColor','none','FontSize',12)
annotation('textbox',[.54 0.63 .4 .2],'String',['Avg. Position (absolute) = ' num2str(meanPositionAbsolute)],'EdgeColor','none','FontSize',12)
set(gca, 'FontSize', 12)

figure
%Error vs regular change in x1
plot(1000*pos_vec(4,:),min_error_vec(1,:),'o')
title(['Variation in error vs change in x1 position (regular) for the first ' num2str(numTrials) ' unflagged trials'])
xlabel('Variation in position (mm)')
ylabel('Normalized Velocity error')
annotation('textbox',[.54 0.7 .4 .2],'String',['Avg. Actual error = ' num2str(meanRegularError)],'EdgeColor','none','FontSize',12)
annotation('textbox',[.54 0.665 .4 .2],'String',['Avg. Optimized error = ' num2str(meanOptimizedError)],'EdgeColor','none','FontSize',12)
annotation('textbox',[.54 0.63 .4 .2],'String',['Avg. Position (regular) = ' num2str(meanPositionRegular)],'EdgeColor','none','FontSize',12)
set(gca, 'FontSize', 12)

figure
%Error vs %change in x1
plot(100*(1-percentage_vec(1,:)),min_error_vec(1,:),'o')
title(['Variation in error vs change (percentage) in x1 position for the first ' num2str(numTrials) ' unflagged trials'])
xlabel('Variation in position (%)')
ylabel('Normalized Velocity error')
annotation('textbox',[.57 0.7 .4 .2],'String',['Avg. Actual error = ' num2str(meanRegularError)],'EdgeColor','none','FontSize',12)
annotation('textbox',[.57 0.665 .4 .2],'String',['Avg. Optimized error = ' num2str(meanOptimizedError)],'EdgeColor','none','FontSize',12)
set(gca, 'FontSize', 12)

% figure
% h = histogram(1000*pos_vec(3,:),30);
% hold on
% plot(conv(h.BinEdges, [0.5 0.5], 'valid'), h.BinCounts)

% Step 7: epsilon/mu analysis
M = zeros(numTrials,3);
for g = 1:numTrials
 M(g,:) = cell2mat(C2{g,1});   
end
AvgEpMu = zeros(2,2);
AvgEpMu(1,1) = mean(M(:,1));
AvgEpMu(2,1) = std(M(:,1));
AvgEpMu(1,2) = mean(M(:,2));
AvgEpMu(2,2) = std(M(:,2));
disp(AvgEpMu)

%Step 8: correlation analysis
figure
%Tangenltial impulse plot
plot(1000*pos_vec(5,:).*impulses(1,:),impulses(1,:),'o')
title({'Optimized tangential impulse vs moment due to tangent', ['The first ' num2str(numTrials) ' unflagged trials'], 'Wang model'})
xlabel('Moment due to tangent')
ylabel('Tangential Impulse')
hold on
plot([xlim],[0 0],'LineWidth',2,'Color','k')
hold on
plot([0 0],[ylim],'LineWidth',2,'Color','k')
set(gca, 'FontSize', 12)

figure
%Normal impulse plot
plot(1000*pos_vec(2,:).*impulses(2,:),impulses(2,:),'o')
title({'Optimized normal impulse vs moment due to normal', ['The first ' num2str(numTrials) ' unflagged trials'], 'Wang model'})
xlabel('Moment due to normal')
ylabel('Normal Impulse')
set(gca, 'FontSize', 12)

figure
%Angle plot
plot(mod(pos_vec(8,:)*180/pi,90),min_error_vec(1,:),'o')
title({'Error vs pre-impact angle', ['The first ' num2str(numTrials) ' unflagged trials'], 'Wang model'})
xlabel('Pre-impact angle')
ylabel('Normalized velocity error')
hold on
plot([90 90],[ylim],'--','LineWidth',2,'Color','k')
hold on
plot([45 45],[ylim],'--','LineWidth',2,'Color','k')
hold on
plot([0 0],[ylim],'--','LineWidth',2,'Color','k')
set(gca, 'FontSize', 12)

% plot(pos_vec(8,:),pos_vec(4,:),'o')
% x = [-8:0.1:8];
% y = sin(2*x)/100;
% hold on
% plot(x,y)

% plot of tangential impulse with the area lines
plot(1000*pos_vec(5,:).*impulses(1,:),impulses(1,:),'.')
title({'Optimized tangential impulse vs moment due to tangent', ['The first ' num2str(numTrials) ' unflagged trials'], 'Wang model'})
xlabel('Moment due to tangent')
ylabel('Tangential Impulse')
hold on
plot([xlim],[0 0],'LineWidth',2,'Color','k')
hold on
plot([0 0],[ylim],'LineWidth',2,'Color','k')
set(gca, 'FontSize', 12)
hold on
x = -40:0.1:30;
%slope1 = abs(min(impulses(1,:))*min(pos_vec(5,:)));
slope1 = 0.029;
y = slope1*x;
plot(x,y,'-','LineWidth',1.5,'Color',[0.4940, 0.1840, 0.5560]);
hold on
%slope2 = abs(min(impulses(1,:))*max(pos_vec(5,:)));
slope2 = 0.0395;
y2 = slope2*x;
plot(x,y2,'-','LineWidth',1.5,'Color',[0, 0.5, 0]);
ylim([-1.2 0.8])
annotation('textbox',[.13 0.68 .4 .2],'String',['Violet Slope = ' num2str(slope1)],'EdgeColor','none','FontSize',12)
annotation('textbox',[.13 0.645 .4 .2],'String',['Green Slope = ' num2str(slope2)],'EdgeColor','none','FontSize',12)

