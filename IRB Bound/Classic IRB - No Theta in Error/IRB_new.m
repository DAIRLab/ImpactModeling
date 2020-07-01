%% Newer IRB

load('ellipse_uniform.mat');

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037;
I1 = m1 * (a0^2 + b0^2) / 4;
Mass = [m1, 0, 0;
    0, m1, 0;
    0, 0, I1];

Tlength = 1000;
errorVec = zeros (1,Tlength);

count  = 0; 

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
        pause(0.5);
        count = count  + 1;
        noWidth(count) = abs(post(3) - predicted(3));
        
        thetaDiff(1, count) = post(3) - predicted(3);
        thetaDiff(2, count) = pre(3) - post(3);
        thetaDiff(3, count) = pre(3) - predicted(3);
        
        
    end
end

avErr =  mean(errorVec);

%% 
plot(thetaDiff(2,:), thetaDiff(1,:), '.');
xlabel("Observed \omega^- - Observed \omega^+")
ylabel("Obsereved \omega^+ - Predicted \omega^+");
title("Running IRB Optimization on Just X_{dot} and Y_{dot}")
%legend("Trial and \omega different Pairs", "Mean of \omega Difference = " + mean(thetaDiff),'Location', 'North')

figure()
plot(thetaDiff(3,:), thetaDiff(1,:), '.');
xlabel("Observed \omega^- - Predicted \omega^+")
ylabel("Observed \omega^+ - Predicted \omega^+");
title("Running IRB Optimization on Just X_{dot} and Y_{dot}")