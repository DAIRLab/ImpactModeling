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
errorP_old = zeros(3,numTrials);

% Step 3: run the for-loop 

for tr = 1:numTrials
    trials = nonflagged(tr);

        p1 = 1;
        p2 = 1;
        %loop over all values
         [P_old,error_old] = IRB_NEW_torque_old(bounce_array,trials,p1,p2);
         [P,w,error] = IRB_NEW_torque(bounce_array,trials,p1,p2);
         errorP(1:3,tr) = P';
         errorP(4,tr) = error;
         errorP(5,tr) = w;
         errorP_old(1:2,tr) = P_old';
         errorP_old(3,tr) = error_old;

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

%Plot torque vs normal impulse*width
plot(errorP(3,:),errorP(5,:).*errorP(2,:),'o')
title({'Torque vs Normal Moment', ['The first ' num2str(numTrials) ' unflagged trials'], 'IRB model'})
xlabel('Normal Moment: Normal P * width (Nms)')
ylabel('Torque (Nms)')
set(gca, 'FontSize', 12)
p = polyfit(errorP(3,:),errorP(5,:).*errorP(2,:),1)
x = -0.003:0.00001:0.005;
y = p(1)*x + p(2);
hold on
plot(x,y,'-','LineWidth',1.5,'Color','r')
plot([xlim],[0 0],'-','LineWidth',1.1,'Color','k')
plot([0 0],[ylim],'-','LineWidth',1.1,'Color','k')
legend('Data Points',['y = ' num2str(p(1)) 'x - ' num2str(abs(p(2)))])

%Comparing torque vs Tangential Moment
tiledlayout(1,2)
nexttile
plot(errorP(1,:).*properties(1,:),errorP(3,:),'.')
title({'Torque vs Tangential Moment', ['The first ' num2str(numTrials) ' unflagged trials'], 'IRB Torque model'})
xlabel('Tangential Moment (M_{IRB}) [Nms]')
ylabel('Torque (M_{W}) [Nms] ')
set(gca, 'FontSize', 12)
p = polyfit(errorP(1,:).*properties(1,:),errorP(3,:),1)
x = -0.002:0.00001:0.0025;
y = p(1)*x + p(2);
hold on
plot(x,y,'-','LineWidth',1.5,'Color','r')
plot([xlim],[0 0],'-','LineWidth',1.1,'Color','k')
plot([0 0],[ylim],'-','LineWidth',1.1,'Color','k')
legend('Data Points',['y = ' num2str(p(1)) 'x + ' num2str(abs(p(2)))])
xlim([-0.002 0.0025])
f1 = polyval(p,errorP(1,:).*properties(1,:));
mean_error1 = mean(f1-errorP(3,:));

nexttile
plot(errorP_old(1,:).*properties(1,:),errorP_old(2,:).*properties(2,:),'.')
title({'Normal Moment vs Tangential Moment', ['The first ' num2str(numTrials) ' unflagged trials'], 'IRB Regular model'})
xlabel('Tangential Moment (M_{IRB}) [Nms]')
ylabel('Normal Moment (M_{W}) [Nms]')
set(gca, 'FontSize', 12)
p2 = polyfit(errorP_old(1,:).*properties(1,:),errorP_old(2,:).*properties(2,:),1)
x = -0.002:0.00001:0.002;
y = p2(1)*x + p2(2);
hold on
plot(x,y,'-','LineWidth',1.5,'Color','r')
plot([xlim],[0 0],'-','LineWidth',1.1,'Color','k')
plot([0 0],[ylim],'-','LineWidth',1.1,'Color','k')
legend('Data Points',['y = ' num2str(p2(1)) 'x + ' num2str(abs(p2(2)))])
f2 = polyval(p2,errorP_old(1,:).*properties(1,:));
mean_error2 = mean(f2-errorP_old(2,:).*properties(2,:));



per_change1 = ((errorP_old(2,:) - errorP(2,:))./errorP_old(2,:))*100;
per_change2 = ((errorP_old(1,:) - errorP(1,:))./errorP_old(1,:))*100;

%Mike's plot
figure
plot(errorP(1,:).*properties(1,:) + errorP(2,:).*properties(2,:),errorP(3,:),'.')
title({'Adjusted Moment vs IRB Moment', ['The first ' num2str(numTrials) ' unflagged trials'], 'IRB Torque model'})
xlabel('M_{IRB} = M_{Tangential} + M_{Normal} [Nms]')
ylabel('M_{W} [Nms]')
set(gca, 'FontSize', 12)
p2 = polyfit(errorP(1,:).*properties(1,:) + errorP(2,:).*properties(2,:),errorP(3,:),1)
x = -0.003:0.00001:0.003;
y = p2(1)*x + p2(2);
hold on
plot(x,y,'-','LineWidth',1.5,'Color','r')
plot([xlim],[0 0],'-','LineWidth',1.1,'Color','k')
plot([0 0],[ylim],'-','LineWidth',1.1,'Color','k')
legend('Data Points',['y = ' num2str(p2(1)) 'x + ' num2str(abs(p2(2)))])

%Angle plots
% plot(mod(properties(3,:)*180/pi,90),errorP_old(3,:),'.')
% hold on
% plot(mod(properties(3,:)*180/pi,360),errorP(4,:),'.')
% ylim([0 0.02])
% plot(errorP(1,:).*properties(1,:) + errorP(3,:),errorP(4,:),'.')
% ylim([0 0.02])

%All moments vs error
figure
plot(errorP(1,:).*properties(1,:) + errorP(2,:).*properties(2,:) + errorP(3,:),errorP(4,:),'.')
title({'All Moments vs Error', ['The first ' num2str(numTrials) ' unflagged trials'], 'IRB Torque model'})
xlabel('M_{Tangential} + M_{Normal} + M_{Torque} [Nms]')
ylabel('Velocity Error Square')
set(gca, 'FontSize', 12)
ylim([0 0.0000000000005])
xlim([-0.0015 0.0015])

%All moments vs width
figure
plot(errorP(1,:).*properties(1,:) + errorP(2,:).*properties(2,:) + errorP(3,:),errorP(3,:),'.')
title({'All Moments vs Torque', ['The first ' num2str(numTrials) ' unflagged trials'], 'IRB Torque model'})
xlabel('M_{Tangential} + M_{Normal} + M_{Torque} [Nms]')
ylabel('M_{Torque} [Nm]')
set(gca, 'FontSize', 12)
%ylim([0 0.05])
xlim([-0.0015 0.0015])
hold on
plot([xlim],[0 0],'-','LineWidth',1.1,'Color','k')
plot([0 0],[ylim],'-','LineWidth',1.1,'Color','k')

%[param]=sineFit(errorP(1,:).*properties(1,:) + errorP(2,:).*properties(2,:) + errorP(3,:),errorP(3,:))

%Output: SineParams(1): offset (offs)
%               SineParams(2): amplitude (amp)
%               SineParams(3): frequency (f)
%               SineParams(4): phaseshift (phi)
hold on
%x2 = [-0.0015:0.00001:0.0015];
x2 = -3.14/1000:0.00001:3.14/1000;
%y = param(2)*sin(param(3)*x2 + param(4)) + param(1);
%y = param(1)+param(2)*sin(2*pi*param(3)*x2+param(4));
amp = -0.0016;
freq = 2600;
shift = 0.75;
up = 0.0004;
y = amp*sin(freq*x2 - shift)+up;
plot(x2,y,'r')
annotation('textbox',[.13 0.01 .4 .2],'String',['Sine Equation = ' num2str(amp) '*sin(' num2str(freq) '*x - ' num2str(shift) ') + ' num2str(up)],'EdgeColor','none','FontSize',12)


%IRB moment vs angular error
%angular difference = post - vcalc
figure
tiledlayout(2,1)
nexttile
plot(errorP(1,:).*properties(1,:) + errorP(2,:).*properties(2,:) + errorP(3,:),errorP(4,:),'.')
title({'All Moments vs Error', ['The first ' num2str(numTrials) ' unflagged trials'], 'IRB Torque model'})
xlabel('M_{Tangential} + M_{Normal} + M_{Torque} [Nms]')
ylabel('Normalized Velocity Error Square')
set(gca, 'FontSize', 12)
ylim([0 0.00000031])
xlim([-0.0015 0.0015])

nexttile
plot(errorP(1,:).*properties(1,:) + errorP(2,:).*properties(2,:) + errorP(3,:),errorP(4,:),'.')
title('Zoomed In')
xlabel('M_{Tangential} + M_{Normal} + M_{Torque} [Nms]')
ylabel('Normalized Velocity (by 100) Error Square')
set(gca, 'FontSize', 12)
ylim([0 0.0000002])
xlim([-0.0015 0.0015])

%angle plot
figure
tiledlayout(3,1)
nexttile
plot(properties(3,:).*(180/pi),errorP(4,:),'.')
title({'Pre-impact Angle vs Error', ['The first ' num2str(numTrials) ' unflagged trials'], 'IRB Torque model'})
xlabel('Pre-impact angle (degrees)')
ylabel('Velocity Error Square')
set(gca, 'FontSize', 12)
ylim([0 5e-13])
hold on
plot([xlim],[0.8e-13 0.8e-13],'-','LineWidth',1.1,'Color','r')
hold on
plot([xlim],[2e-13 2e-13],'-','LineWidth',1.1,'Color','r')
hold on
plot([min(xlim) min(xlim)],[0.8e-13 2e-13],'-','LineWidth',1.1,'Color','r')
hold on
plot([max(xlim) max(xlim)],[0.8e-13 2e-13],'-','LineWidth',1.1,'Color','r')

nexttile
plot(properties(3,:).*(180/pi),errorP(4,:),'.')
title('Zoomed In')
xlabel('Pre-impact angle (degrees)')
ylabel('Velocity Error Square')
set(gca, 'FontSize', 12)
ylim([0.8e-13 2e-13])
hold on
plot([xlim],[min(ylim) min(ylim)],'-','LineWidth',1.1,'Color','r')
hold on
plot([xlim],[max(ylim) max(ylim)],'-','LineWidth',1.1,'Color','r')
hold on
plot([min(xlim) min(xlim)],[0.8e-13 2e-13],'-','LineWidth',1.1,'Color','r')
hold on
plot([max(xlim) max(xlim)],[0.8e-13 2e-13],'-','LineWidth',1.1,'Color','r')
hold on
plot([xlim],[0.8e-13 0.8e-13],'-','LineWidth',1.1,'Color','g')
hold on
plot([xlim],[0.9e-13 0.9e-13],'-','LineWidth',1.1,'Color','g')
hold on
plot([min(xlim) min(xlim)],[0.8e-13 0.9e-13],'-','LineWidth',1.1,'Color','g')
hold on
plot([max(xlim) max(xlim)],[0.8e-13 0.9e-13],'-','LineWidth',1.1,'Color','g')

nexttile
plot(properties(3,:).*(180/pi),errorP(4,:),'.')
title('Even more Zoomed In')
xlabel('Pre-impact angle (degrees)')
ylabel('Velocity Error Square')
set(gca, 'FontSize', 12)
ylim([0.8e-13 0.9e-13])
hold on
plot([xlim],[min(ylim) min(ylim)],'-','LineWidth',1.1,'Color','g')
hold on
plot([xlim],[max(ylim) max(ylim)],'-','LineWidth',1.1,'Color','g')
hold on
plot([min(xlim) min(xlim)],[0.8e-13 0.9e-13],'-','LineWidth',1.1,'Color','g')
hold on
plot([max(xlim) max(xlim)],[0.8e-13 0.9e-13],'-','LineWidth',1.1,'Color','g')
