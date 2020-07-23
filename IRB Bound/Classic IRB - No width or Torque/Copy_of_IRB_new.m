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

Tlength = 92;
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

    end
end

avErr =  mean(errorVec, 2);
disp(avErr);

%%
figure()
plot(useful(1,:), errorVec, '.')
ylabel("Normalized l2 Norm Velocity Error");
xlabel("Pre-Impact Angle");
title("IRB No Torque, Square Data");

figure()
plot(1:count, errorVec, '.');
hold on
ylabel("Normalized l2 Norm Velocity Error");
xlabel("Trial");
title("IRB No Torque, Square Data");


