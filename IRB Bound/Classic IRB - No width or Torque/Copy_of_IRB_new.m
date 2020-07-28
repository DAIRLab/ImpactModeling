%% Newer IRB with Square Data
clear;
load('squareDataPhil.mat');

% Set up Constants
stepSize = 0.01;
sl = 0.06; %side length of square from data README
m1 = 0.0485; %WHAT IS THE MASS????? Right now we have an area ratio 
                                    %approximation from ellipse

I1 = m1 * sl^2 / 6; %Moment of Inertia

Mass = [m1, 0, 0;
        0, m1, 0;
        0, 0, I1]; %generlaized Mass matrix

Tlength = 387;
errorVec = zeros (1,Tlength);

count  = 0; 

for i = 1:Tlength

    %Select trial
    trial = i;

    if true%sum(bounce_array(trial).flags) < 1
        % Finding pre and post impact velocities / states
        pre = squareDataPhil(trial).states(4:6)';
        post = squareDataPhil(trial).states(10:12)';

        d = (squareDataPhil(trial).d);   %tangential
        n = (squareDataPhil(trial).n);   %normal

        J = [d;n]; %Jacobian

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

        error = findError(P, Mass, J, pre, post); %final error
        predicted = pre + inv(Mass) * J' * P'; %predicted post impact state

        count = count  + 1;
        errorVec(count) = error;

%         disp("Observed:")
%         disp(post)
%         disp("Predicted:")
%         disp(predicted);

        %vector for keeping track of variables to plot/look for
        %correlations
        useful(1,count) = squareDataPhil(trial).states(3);
        useful(2,count) = abs(post(3) - pre(3));
        useful(3, count) = abs(post(3));


    end
end

avErr =  mean(errorVec, 2);

figure();
plot(useful(2,:), errorVec, '.')
xlabel("Change in Angular Velocity");
ylabel("Normalized l2 Norm Velocity Error");

%%
for j = 1:count
    wrappedAngle(j) = wrapTo180(rad2deg(useful(1,j)));
end
figure()
hold on
plot(1:count, wrappedAngle, '.')
plot(linspace(0,count), zeros(1, 100), 'r-')
plot(linspace(0,count), 90* ones(1, 100), 'r-')
plot(linspace(0,count), -90 * ones(1,100), 'r-')

ylabel("Wrapped Pre-Impact Angle");
xlabel("Trial Number");
title("Square Data");
ylim([-180, 180]);
xlim([0, count]);
%%
figure()
plot(wrappedAngle, errorVec, '.');
hold on
ylabel("Normalized l2 Norm Velocity Error");
xlabel("Wrapped Pre-Impact Angle");
%title("IRB No Width, Square Data");

figure()
plot(useful(1,:), errorVec, '.')
ylabel("Normalized l2 Norm Velocity Error");
xlabel("Pre-Impact Angle");
title("IRB No Torque, Square Data");
