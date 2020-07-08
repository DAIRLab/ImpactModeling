%% Newer IRB
clear; 

load('ellipse_uniform.mat');

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
Mass = [m1, 0, 0;
        0, m1, 0; 
        0, 0, I1]; 
    

widths = [0.1];%linspace(0, 0.01, 1);
numTrials = 2000; %Number of Trials
ran = randi([1 2000], 1, numTrials);

count = 0; %counter variable
bestWidth = [];
for a = 1:length(widths) 
    %choose the current max allowable width
    maxWidth = widths(a);
    %reset the error vector matrix
    errorVector = [];
    
    for i = 1:length(ran)
        trial = i; %ran(i);
        if sum(bounce_array(trial).flags) < 1
            if (count == 167) 
                disp(i);
            end
            % Finding pre and post impact velocities / states
            pre = bounce_array(trial).states(4:6)';
            post = bounce_array(trial).states(10:12)';

            d = (bounce_array(trial).d);   %tangential
            n = (bounce_array(trial).n);   %normal

            J = [d;n; 0, 0, 1]; %Jacobian

            fun = @(P)(findError(P, Mass, J, pre, post));
            nonlcon = @(P)(constraint(P, Mass, J, maxWidth, pre, post));

            P0 = [0 0 0];
            A = []; % No other constraints
            b = [];
            Aeq = [];
            beq = [];
            lb = [];
            ub = [];
            options = optimoptions('fmincon','FiniteDifferenceType','central', ...
                    'StepTolerance',1e-10, 'Display','off');
            P = fmincon(fun, P0, A, b, Aeq, beq, lb, ub, nonlcon, options);
            %iterator for non flagged trials 
            count = count + 1;
            %save some data for correlation plots
            bestWidth(1, count) = P(1);
            bestWidth(2,count) = P(2);
            bestWidth(3,count) = P(3);
%             bestWidth(4, count) = 1/2 * pre' * Mass * pre;
%             bestWidth(5, count) = 1/2 * post' * Mass * post;
%             bestWidth(6, count) = bounce_array(trial).states(3);
%             bestWidth(7, count) = bounce_array(trial).states(6);
            %ellipse_visual(pre(1), pre(2), pre(3), 'b');
            %add the trial's error to error vector
            error = findError(P, Mass, J, pre, post);
            errorVec(count) = error;
%             disp("Trial: " + trial);
%             disp("Pre Impact Angle: " + (rem(bounce_array(trial).states(3), pi)*180)/pi);
%             disp("Post Impact Omega Dot: " + post(3));
            predicted = pre + inv(Mass) * J' * [P(1:2)'; P(3) * P(2)];
            
%             disp("Predicted Post Impact Omega Dot: " + predicted(3));
%             disp("Optimal Width: " + P(3));

              yesWidth(1,count) = abs(post(3) - predicted(3));
              yesWidth(2, count) = d(3) * P(1) + n(3) * P(2);
%             
%             checkWidth(1, count) = abs(post(3)) - abs(predicted(3));
%             checkWidth(2, count) = error;
%             posaTest = J' * [P(1:2)'; P(3) * P(2)];
%             checkWidth(3, count) = posaTest(3);
            
         end

    end
    
    %find the average error across all trials for a certain max width
    avgError(a) = mean(errorVector);
end

% disp("Tangential Impulse: " + P(1) + " [N*s]")
% disp("Normal Impluse: " + P(2) + " [N*s]")

%% Look at error for IRB with w/ without a max constraint, should be nearly 0
errorVec(1420) = 0;
errorVec(310) = 0;
errorVec(234) = 0;

plot(1:count, errorVec(:), '.')
xlabel("Unflagged Trial #")
ylabel("(Normalized Error)^2")
title("IRB w/ Width but no Max Width Constraint")



%%




figure()
plot(yesWidth(1,:), yesWidth(2,:), '.');
xlabel("|\omega_{observed} - \omega_{predicted}| ");
ylabel("IRB Moment");
title("Angular Velocity Error Vs. IRB Moment");
%%
figure()
plot(bestWidth(1, :), bestWidth(3, :), '.');
ylabel("Optimal Width");
xlabel("Tangential Impulse");
xlim([-0.03, 0.03])
ylim([-0.01, 0.01])
title("Not Squaring Error");
%%
clumps = 0;
for i = 1:length(bestWidth)
    if(abs(bestWidth(3,i)) < 0.0001)
        clumps = clumps  + 1;
        %disp(i);
    end
    if abs((bestWidth(1,i) - 0.01188)) < 0.0001
        disp(i)
    end
end

disp(clumps)

% 85 squaring and 59 not squaring
%% Post Process Data
%avgW = max(bestWidth, [], 1);
%figure();
%plot(widths, avgError)
%title("Maximum Allowable Width vs. Error Plot")
%xlabel("Maximum Allowable Width [m]")
%ylabel("Average Normalized Error Across All Trials")

% figure();
% hold on
% for j = 1:100
%     plot(bounce_array(j).states(4), avgW(j), 'r*');
%     disp(j)
% end 
% 
% ylabel("Avg Optimal Width [m]")
% xlabel("Pre-Impact X Velocity [m]")
% title("Optimal Width and X Velocity Correlation")
% figure();
% hold on
% for j = 1:100
%     plot(bounce_array(j).states(5), avgW(j), 'b*');
%     disp(j)
% end 
% title("Optimal Width and Y Velocity Correlation")
% ylabel("Avg Optimal Width [m]")
% xlabel("Pre-Impact Y Velocity [m]")
% 
figure();
hold on
for j = 1:count
    theta(j) = bestWidth(6,j);
    if abs(theta(j))> pi
        theta(j) = (rem(theta(j),pi)*180)/pi;
    else
        theta(j) = theta(j)*180/pi;
    end
    %plot(theta, bestWidth(3, j), 'b.');
end 
plot(theta, bestWidth(3, :), 'b.');
plot(linspace(-180,180), (sind(2*linspace(-180,180))+0.3)/100, 'r')
%title("Optimal Width and Impact Angle Correlation")
ylabel("Avg Optimal Width [m]")
xlabel("Pre-Impact Angle [degrees]")
xlim([-180, 180]);
ylim([-0.02, 0.02]);
legend("Optimal Width Angle Pairs", "(Sin(2\theta) + 0.3)/100", "Location", "Southwest")
%%
figure();
hold on
plot(bestWidth(2,:), bestWidth(3,:), '.');
plot(ones(1, 100) * 0.1086, linspace(-0.01, 0.01));
legend('Optimal Width Normal Impulse Pairs', 'Minimum Impulse For Nonzero Optimal Width', ...
    'Location', 'Northwest');
xlabel("Normal Impulse")
ylabel("Optimal Width [m]")
title("Optimal Width and Normal Impulse Correlation")

var1 = [bestWidth(2,:)', bestWidth(3,:)'];
Q = corrcoef(var1)
var2 = [bestWidth(1,:)', bestWidth(3,:)'];
R = corrcoef(var2)

figure()
plot(bestWidth(1,:), bestWidth(3,:), '.');

xlabel("Tangential Impulse[N*s]")
ylabel("Optimal Width [m]")
%title("Optimal Width and Tangential Impulse Correlation")

figure()
plot(bestWidth(4,:), bestWidth(3,:), '.');

xlabel("Pre Impact Energy [J]")
ylabel("Optimal Width [m]")
%title("Optimal Width and Tangential Impulse Correlation")

figure()
plot(bestWidth(5,:), bestWidth(3,:), '.');

xlabel("Post Impact Energy [J]")
ylabel("Optimal Width [m]")
%title("Optimal Width and Tangential Impulse Correlation")

figure()
plot(bestWidth(7,:), bestWidth(3,:), '.');

xlabel("Pre Impact \dot\theta [rad/s]")
ylabel("Optimal Width [m]")
%title("Optimal Width and Tangential Impulse Correlation")



% %% 
% nonZ = [];
% for i = 1:length(bestWidth(1,:))
%    if(abs(bestWidth(3,i)) > 0.0005)
%        nonZ = [nonZ, bestWidth(2,i)];
%    end
% end
% min(nonZ)
