%plot various graphs to show correlations between error, impulse, moment

load('ellipse_uniform.mat');

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.0364; 
I1 = m1 * (a0^2 + b0^2) / 4;
Mass = [m1, 0, 0;
        0, m1, 0; 
        0, 0, I1]; 
numTrials = 1000;
ran = randi([1 2000], 1, numTrials);
errors = [];
width = [];
P_n = [];
P_t = [];
moment_n = [];
moment_t = [];
moment_net = [];
count = 0; 

for i = 1:numTrials
    z = ran(i);
    if sum(bounce_array(z).flags) == 0 
    % Finding pre and post impact velocities / states
    pre = bounce_array(z).states(4:6)';
    post = bounce_array(z).states(10:12)';
    yPos = bounce_array(z).states(2);

    d = (bounce_array(z).d);   %tangential
    n = (bounce_array(z).n);   %normal

    J = [d; n; 0 0 .5]; %Jacobian
    
    fun = @(P)(findError_torque(P, Mass, J, pre, post));
    nonlcon = @(P)(constraint_torque(P, Mass, J, pre, post));

    P0 = [0 0 0]; %Pt Pn torque
    A = []; % No other constraints
    b = [];
    Aeq = [];
    beq = [];
    lb = [];
    ub = [];
    P = fmincon(fun, P0, A, b, Aeq, beq, lb, ub, nonlcon);  
    
    w = abs(P(3)) / P(2);

    width = [width; w];
    errors = [errors; findError_torque(P, Mass, J, pre, post)];
    P_n = [P_n; P(2)];
    P_t = [P_t; P(1)];
    moment_n = [moment_n; P(3)];
    moment_t = [moment_t; P(1)*yPos]; %this is assuming that the pre-impact y coordinate is the y-coordinate of the COM at impact **WILL NOT give highly accurate results
    moment_net = [moment_net; moment_n(end) + moment_t(end)];
    count = count + 1;
    
    end
end

figure(1)
scatter(1000*width, errors,'filled', 'b')
title("Patch Size vs Error for " + count + " Random Trials")
xlabel('Width [mm]')
ylabel('Normalized Error')

figure(2)
scatter(1000*moment_n, P_n, 'filled', 'b')
hold on
line([0 0], [0.1 0.2], 'Color', 'black');
title("Moment due to Normal Impulse vs Normal Impulse for " + count + " Random Trials")
xlabel('Moment (Normal) [Nmm]')
ylabel('Normal Impulse [N*s]')

figure(3)
scatter(1000*moment_n, P_t, 'filled', 'b')
hold on
line([-5 5],[0 0], 'Color', 'black'); 
hold on
line([0 0], [-.04 0.03], 'Color', 'black');
title("Moment due to Normal Impulse vs Tangential Impulse for " + count + " Random Trials")
xlabel("Moment (Normal) [Nmm]")
ylabel("Tangential Impulse [N*s]")

figure(4)
scatter(1000*moment_t, P_n, 'filled','b')
hold on
line([0 0], [0.1 0.2], 'Color', 'black');
title("Moment due to Tangential Impulse vs Normal Impulse for " + count + " Random Trials")
xlabel("Moment (Tangential) [Nmm]")
ylabel("Normal Impulse [N*s]")

figure(5)
scatter(1000*moment_t, P_t, 'filled', 'b')
hold on
line([-1.5 1],[0 0], 'Color', 'black'); 
hold on
line([0 0], [-.04 0.03], 'Color', 'black');
title("Moment due to Tangential Impulse vs Tangential Impulse for " + count + " Random Trials")
xlabel("Moment (Tangential) [Nmm]")
ylabel("Tangential Impulse [N*s]")

figure(6)
scatter(1000*moment_net, P_n, 'filled', 'b')
hold on
line([0 0], [0.1 0.2], 'Color', 'black');
title("Net Moment vs Normal Impulse for " + count + " Random Trials")
xlabel("Net Moment [Nmm]")
ylabel("Normal Impulse [N*s]")

figure(7)
scatter(1000*moment_net, P_t, 'filled', 'b')
hold on
line([-5 4],[0 0], 'Color', 'black'); 
hold on
line([0 0], [-.04 0.03], 'Color', 'black');
title("Net Moment vs Tangential Impulse for " + count + " Random Trials")
xlabel("Net Moment [Nmm]")
ylabel("Tangential Impulse [N*s]")
