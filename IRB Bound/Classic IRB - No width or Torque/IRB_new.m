%% Newer IRB
clear;
load('ellipse_uniform.mat');

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037;
I1 = m1 * (a0^2 + b0^2) / 4;

spacer =  logspace(-2, 2, 30);%100, 50, 10, 5, 2, 1, 0.5, 0.2];
Tlength = 2000;
%errorVec = zeros (1,Tlength);
for a = 1%:length(spacer)
    I1 = m1 * (a0^2 + b0^2) / 4;
    Mass = [m1, 0, 0;
    0, m1, 0;
    0, 0, I1];
    count  = 0; 
    scale = spacer(a);
    for i = 1:Tlength

        %Select trial
        trial = i;

        if sum(bounce_array(trial).flags) < 1
            % Finding pre and post impact velocities / states
            pre = bounce_array(trial).states(4:6)';
            post = bounce_array(trial).states(10:12)';

            d = (bounce_array(trial).d);   %tangential
            n = (bounce_array(trial).n);   %normal

            J = [d;n]; %Jacobian
            
            fun = @(P)(findError(P, Mass, J, pre, post));
            %fun = @(P)(findErrorScaled(P, Mass, J, pre, post, scale));
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

            error = findError(P, Mass, J, pre, post);

    %         disp("Trial: " + i);
    %         disp("Pre Impact Angle: " + (rem(bounce_array(i).states(3), pi)*180)/pi);
    %         disp("Post Impact Omega Dot: " + post(3));
             predicted = pre + inv(Mass) * J' * P';
    %         disp("Predicted Post Impact Omega Dot: " + predicted(3));
    %         disp("Difference: " + abs(post(3) - predicted(3)));
            %disp("Optimal Width: " + P(3));
            count = count  + 1;
            errorVec(a, count) = error;
            %noWidth(count) = abs(post(3) - predicted(3));

            %checkNoWidth(1, count) = abs(post(3)) - abs(predicted(3));
            %checkNoWidth(2, count) = error;
            yesWidth(1,count) = abs(post(3) - predicted(3));
            yesWidth(2, count) = d(3) * P(1) + n(3) * P(2);
            useful(1, count) = wrapTo180(rad2deg(bounce_array(trial).states(3)));
            useful(2,count) = abs(post(3) - pre(3));
            useful(3, count) = abs(post(3));

        end
    end
end

avErr =  mean(errorVec, 2);
disp(avErr);
figure()
plot(useful(3,:), errorVec, '.')
xlabel("Post Impact Angular Velocity");
ylabel("Normalized l2 Norm Velocity Error");
%%
figure()
plot(useful, errorVec, '.')
xlabel("Wrapped Pre-Impact Angle")
ylabel("Normalized l2 Norm Velocity Error");


%%
figure()
semilogx(spacer, avErr);
hold on
%semilogx(spacer(avErr == min(avErr)), min(avErr), 'r*');
ylabel("Normalized l2 Norm Velocity Error");
xlabel("Scale in Cost Function for \theta_{dot}");
%title("IRB No Torque");


%%
figure()
plot(1:count, errorVec(:), '.')
xlabel("Unflagged Trial #")
ylabel("Error (Using Original Metric)");
title("Optimizing with Scaled Cost Function");


% figure()
% plot(yesWidth(1,:), yesWidth(2,:), '.');
% xlabel("|\omega_{observed} - \omega_{predicted}| ");
% ylabel("IRB Moment");
% title("Angular Velocity Error Vs. IRB Moment");
% 
% 
% figure()
% plot(1:count, errorVec(3,:), '.')
% xlabel("Trial")
% ylabel("Normalized Veloicty Error")
% title("Normal IRB (No Torque)");

%%
figure()
plot(1:2.5:50, avErr')
xlabel("Percent Decrease in Moment of Inertia");
ylabel("Average Error Across Non-Flagged Trials")
%title("Moment of Inertia vs Average Error");
