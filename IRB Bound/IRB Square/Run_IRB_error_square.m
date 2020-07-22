% Run_IRB_error
        
numTrials = 41;
errorP = zeros(5,numTrials);
errorP_old = zeros(3,numTrials);

% Step 3: run the for-loop 

for tr = 1:numTrials

        p1 = 1;
        p2 = 1;
        %loop over all values
         [P_old,error_old] = IRB_NEW_torque_old_square(squareData,tr,p1,p2);
         [P,w,error] = IRB_NEW_torque_square(squareData,tr,p1,p2);
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
    
    
    D = squareData(tr).d;
    N = squareData(tr).n;
    ang = squareData(tr).states(5);

    properties(1,tr) = D(3); 
    properties(2,tr) = N(3);
    properties(3,tr) = ang;

end

% Step 5: typical plots
%Moment(tangent) vs Tangential impulse


IRB_moment = errorP(1,:).*properties(1,:) + errorP(2,:).*properties(2,:);

%IRB vs Torque
plot(errorP(1,:).*properties(1,:) + errorP(2,:).*properties(2,:),errorP(3,:),'o')
p = polyfit(errorP(1,:).*properties(1,:) + errorP(2,:).*properties(2,:),errorP(3,:),1)
x = -0.008:0.00001:0.003;
y = p(1)*x + p(2);
hold on
plot(x,y,'-','LineWidth',1.5,'Color','r')
figure
%All moments vs Error
plot(errorP(1,:).*properties(1,:) + errorP(2,:).*properties(2,:) + errorP(3,:),errorP(4,:),'o')
ylim([0 10e-13])
figure
%Error vs angle
plot(properties(3,:),errorP(4,:),'o')
ylim([0 10e-13])

