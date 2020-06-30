% Run_IRB_error

% Step 1: load in ellipse data 
load('ellipse_uniform.mat');

% Step 2: set-up the number of trials
desiredTrials = 2000;
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
errorP = zeros(5,numTrials);

% Step 3: run the for-loop 

for tr = 1:numTrials
    trials = nonflagged(tr);

        %loop over all values
         %[P,error] = IRB_bound2(bounce_array,trials);
         [P,w,error] = IRB_NEW_torque(bounce_array,trials);
         errorP(1:3,tr) = P';
         errorP(4,tr) = error;
         errorP(5,tr) = w;

    disp(tr)   
end

% Step 4: obtain trial properties
properties = zeros(3,numTrials);

for tr = 1:numTrials
    trials = nonflagged(tr);
    D = bounce_array(trials).d; %y
    N = bounce_array(trials).n; %x
    ang = bounce_array(trials).states(3);

    properties(1,tr) = D(3); 
    properties(2,tr) = N(3);
    properties(3,tr) = ang;

end

% Step 5: typical plots
%Moment(tangent) vs Tangential impulse
plot(errorP(1,:).*properties(1,:),errorP(1,:),'.')
%Moment(normal) vs Normal impulse
figure
plot(errorP(2,:).*properties(2,:),errorP(2,:),'.')
%Angle vs error
figure
plot(mod(properties(3,:)*180/pi,360),errorP(4,:),'.')
ylim([0 0.1])
%Total moment vs shift
IRB_moment = errorP(1,:).*properties(1,:) + errorP(2,:).*properties(2,:);
%total_moment = IRB_moment + errorP(3,:);
total_moment = errorP(1,:).*properties(1,:) + errorP(3,:);
figure
plot(total_moment,errorP(5,:).*errorP(2,:),'.')
title({'Net Moment vs Optimal Patch Size', ['The first ' num2str(numTrials) ' unflagged trials'], 'IRB model'})
xlabel('Tangential moment + Torque (Nms)')
ylabel('Width*Normal Impulse (Nms)')
set(gca, 'FontSize', 12)
p = polyfit(total_moment,errorP(5,:).*errorP(2,:),1)
x = -0.002:0.00001:0.0015;
y = p(1)*x + p(2);
hold on
plot(x,y,'-','LineWidth',1.5,'Color','r')
plot([xlim],[0 0],'-','LineWidth',1.1,'Color','k')
plot([0 0],[ylim],'-','LineWidth',1.1,'Color','k')
legend('Data Points',['y = ' num2str(p(1)) 'x - ' num2str(abs(p(2)))])

