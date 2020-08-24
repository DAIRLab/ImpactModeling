%% Explore Sensitivity of Jacobian

%Run IRB without torque on trials with high error and allow for a certain
%percentage of angle adjustment, choose the new angle to be that with the
%best prediction due to IRB and note what change was made (percentage wise)


clear;
load('ellipse_uniform.mat');
load('highError.mat');


% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037;
Rg = sqrt(a0^2 + b0^2)*0.5;
I1 = m1 * (a0^2 + b0^2) / 4;

Mass = [m1, 0, 0;
0, m1, 0;
0, 0, I1];

count = 0;
Tlength = length(highError);
percentages = -0.01:0.001:0.01;

for i = 1:Tlength

        %Select trial
        trial = highError(i, 1);
        
        if sum(bounce_array(trial).flags) < 1
            minError = 100; 
           % Finding pre and post impact velocities / states
            pre = bounce_array(trial).states(4:6)';
            post = bounce_array(trial).states(10:12)';

            for j = 1:length(percentages)
                
                %Find jacobian using Andy's method and updated angle (from
                %percentages)
                angle = bounce_array(trial).states(3) + ... 
                    percentages(j) * bounce_array(trial).states(3);
                
                J = getJacobian(bounce_array(trial).states(1), ...
                bounce_array(trial).states(2), angle);

                
                %Run IRB model with certian Jacobian
                fun = @(P)(findError(P, Mass, J, pre, post));
                nonlcon = @(P)(constraint(P, Mass, J, pre, post));

                P0 = [0 0];
                A = []; % No other constraints
                b = [];
                Aeq = [];
                beq = [];
                lb = [];
                ub = [];
                options = optimoptions('fmincon','FiniteDifferenceType','central', ...
                            'StepTolerance',1e-10, 'Display','off');

                P = fmincon(fun, P0, A, b, Aeq, beq, lb, ub, nonlcon, options);

                normError = findError(P, Mass, J, pre, post);
                error = normError*norm([post(1:2);Rg*post(3)]);
                
                if (error < minError)
                    bestPercentage = percentages(j);
                    minError = error;
                end
                
                
            end

        count = count  + 1;
        out(count, :) = [bestPercentage, highError(i, 2), minError];

        end
end
%% 
figure();
hold on
plot(1:count, out(:, 2),'.');
plot(1:count, out(:, 3), '.');
xlabel("Trial")
ylabel("Scaled l2 Norm Velocity Error");
legend("Original Error (Classic IRB)", "New Error with up to +/-1% Change in Angle")
