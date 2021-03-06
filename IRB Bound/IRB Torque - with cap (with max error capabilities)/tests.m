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
    
    for i = 2 %1:length(ran)
        trial = i; %ran(i);
        if sum(bounce_array(trial).flags) < 1
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
            
            Copy_of_constraint(P, Mass, J, maxWidth, pre, post);
            %save some data for correlation plots
%             bestWidth(1, count) = P(1);
%             bestWidth(2,count) = P(2);
%             bestWidth(3,count) = P(3);
%             bestWidth(4, count) = 1/2 * pre' * Mass * pre;
%             bestWidth(5, count) = 1/2 * post' * Mass * post;
%             bestWidth(6, count) = bounce_array(trial).states(3);
%             bestWidth(7, count) = bounce_array(trial).states(6);
            %ellipse_visual(pre(1), pre(2), pre(3), 'b');
            %add the trial's error to error vector
            error = findError(P, Mass, J, pre, post);
            errorVec(count) = error;

            predicted = pre + inv(Mass) * J' * [P(1:2)'; P(3) * P(2)];
            
         end

    end
    
    %find the average error across all trials for a certain max width
    avgError(a) = mean(errorVector);
end
figure()
hold on 
Plot_energy_Ellipse(2);
plot(P(1), P(2), 'r*');


% disp("Tangential Impulse: " + P(1) + " [N*s]")
% disp("Normal Impluse: " + P(2) + " [N*s]")
