%% Newer IRB
clear;
load('ellipse_uniform.mat');

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
Rg = sqrt(a0^2 + b0^2)*0.5;
m1 = 0.037;
I1 = m1 * (a0^2 + b0^2) / 4;
highError = [];

Tlength = 2000;
%errorVec = zeros (1,Tlength);
Mass = [m1, 0, 0;
0, m1, 0;
0, 0, I1];
count  = 0; 
ct = 0;
nima = 0; 
juniors = 0;
ken = 0;
kej = 0;

for i = 1:Tlength

    %Select trial
    trial = i;

    if sum(bounce_array(trial).flags) < 1
        count = count  + 1;
        % Finding pre and post impact velocities / states
        pre = bounce_array(trial).states(4:6)';
        post = bounce_array(trial).states(10:12)';

        d = (bounce_array(trial).d);   %tangential
        n = (bounce_array(trial).n);   %normal
        
        P0 = [0 0];
        A = []; % No other constraints
        b = [];
        Aeq = [];
        beq = [];
        lb = [];
        ub = [];
        options = optimoptions('fmincon','FiniteDifferenceType','central', ...
                    'StepTolerance',1e-10, 'Display','off');

        %Choose Jacobian
        for j = 1:2
            if j == 1
                J = [d;n]; %Jacobian
                out(count, 1:2) = J(:,3)';
            else 
                J = getJacobian(bounce_array(trial).states(1:3));
                out(count, 3:4) = J(:,3)';
            end
            
            fun = @(P)(findError(P, Mass, J, pre, post));
            %fun = @(P)(findErrorScaled(P, Mass, J, pre, post, scale));
            nonlcon = @(P)(constraint(P, Mass, J, pre, post));

            P_out(j,:) = fmincon(fun, P0, A, b, Aeq, beq, lb, ub, nonlcon, options);
            
            error(j) = findError(P_out(j, :), Mass, J, pre, post);
        end

        errorVec(count) = min(error);%*norm([post(1:2);Rg*post(3)]);
        if find(error == min(error)) == 1 
            nima = nima + 1;
            ken = ken + 1/2 * pre' * Mass * pre;
            out(count, 5) = 1;
            predicted = pre + inv(Mass) * [d;n]' * P_out(1,:)';
            out(count, 6) = sign(post(1));
            out(count, 7) = sign(predicted(1));
            out(count, 8) = abs(post(3) - predicted(3));
        else 
            juniors  = juniors + 1;
            kej = kej + 1/2 * pre' * Mass * pre;
            out(count, 5) = 2;
            J = getJacobian(bounce_array(trial).states(1:3));
            predicted = pre + inv(Mass) * J' * P_out(2,:)';
            out(count, 6) = sign(post(1));
            out(count, 7) = sign(predicted(1));
            out(count, 8) = abs(post(3) - predicted(3));
        end
        
        useful(1, count) = wrapTo180(rad2deg(bounce_array(trial).states(3)));
        useful(2,count) = abs(post(3) - pre(3));
        useful(3, count) = abs(post(3));
        useful(4, count) = norm(post);
        useful(5, count) = norm(pre);
    end
end

%%
disp("Nima:")
disp(nima)
disp("Avg Pre Impact Energy:")
disp(ken/nima);
disp("Juniors")
disp(juniors)
disp("Avg Pre Impact Energy:")
disp(kej/juniors);

text = "Using Jacobian w/ Lowest Error";

figure
plot(useful(2,:), errorVec, '.')
xlabel("Change In Angular Velocity")
ylabel("Scaled l2 Norm Velocity Error");
title(text)

figure
plot(useful(1,:), errorVec, '.')
xlabel("Pre Impact Wrapped Angle")
ylabel("Scaled l2 Norm Velocity Error");
title(text)

figure
hold on
plot(1:count, errorVec, '.')
plot(linspace(1,count), ones(1,100)*median(errorVec))
plot(linspace(1,count), ones(1,100)*mean(errorVec))
xlabel("Trial Number")
ylabel("Scaled l2 Norm Velocity Error");
t1 = ["Median Error: " + num2str(median(errorVec))];
t2 = ["Mean Error: " + num2str(mean(errorVec))];
legend("Individual Trial Error", t1, t2);
%title(text)

%%
figure()
plot(useful, errorVec, '.')
xlabel("Wrapped Pre-Impact Angle")
ylabel("l2 Norm Velocity Error");


