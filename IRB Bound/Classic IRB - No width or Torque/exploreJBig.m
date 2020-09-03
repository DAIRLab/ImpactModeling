load('ellipse_uniform.mat');

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
Rg = sqrt(a0^2 + b0^2)*0.5;
m1 = 0.037;
I1 = m1 * (a0^2 + b0^2) / 4;
highError = [];

%errorVec = zeros (1,Tlength);
Mass = [m1, 0, 0;
0, m1, 0;
0, 0, I1];
count  = 0; 
Tlength = 2000;

for j = 1:Tlength
    trial = j;
    errorVec = [];
    if sum(bounce_array(trial).flags) < 1
        count = count + 1;
        % Finding pre and post impact velocities / states
        pre = bounce_array(trial).states(4:6)';
        post = bounce_array(trial).states(10:12)';

        d = (bounce_array(trial).d);   %tangential
        n = (bounce_array(trial).n);   %normal

        J_nima = [d;n]; %Jacobian
        J_new = getJacobian(bounce_array(trial).states(1:3));

        nice = linspace(min([J_new(2,3), J_nima(2,3)]), max([J_new(2,3), J_nima(2,3)]), 20);
        for i = 1:length(nice)
            J = J_nima;
            J(2,3) = nice(i);

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

            errorVec(i) = error; 
        end
        bigErrorVec(count, :) = [min(errorVec), find(errorVec == min(errorVec))];
    end
end
%%
figure
hold on
plot(1:count, bigErrorVec(:,1), '.')
plot(linspace(1,count), ones(1,100).*median(bigErrorVec(:,1)))
plot(linspace(1,count), ones(1,100).*mean(bigErrorVec(:,1)))

xlabel("Trial");
ylabel("Error");
title("Using Optimal J(2,3)")
t1 = ["Median Error: " + num2str(median(bigErrorVec(:,1)))];
t2 = ["Mean Error: " + num2str(mean(bigErrorVec(:,1)))];
legend("Individual Trial Error", t1, t2);