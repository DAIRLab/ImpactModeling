load('ellipse_uniform.mat');

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
Rg = sqrt(a0^2 + b0^2)*0.5;
m1 = 0.037;
I1 = m1 * (a0^2 + b0^2) / 4;

Tlength = 2000;
%errorVec = zeros (1,Tlength);
Mass = [m1, 0, 0;
0, m1, 0;
0, 0, I1];
count  = 0; 
for j = 1:length(highError)
    trial = highError(j);
    errorVec = [];
    % Finding pre and post impact velocities / states
    pre = bounce_array(trial).states(4:6)';
    post = bounce_array(trial).states(10:12)';

    d = (bounce_array(trial).d);   %tangential
    n = (bounce_array(trial).n);   %normal

    J_nima = [d;n]; %Jacobian
    J_new = getJacobian(bounce_array(trial).states(1:3));

    nice = linspace(min([J_new(2,3), J_nima(2,3)]), max([J_new(2,3), J_nima(2,3)]), 50);
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
    clf;
    hold on
    plot(nice, errorVec, '.')
    plot(ones(10)*J_nima(2,3), linspace(min(errorVec), max(errorVec),10), 'r');
    plot(ones(10)*J_new(2,3), linspace(min(errorVec), max(errorVec), 10), 'g');
    xlabel("J(2,3)");
    ylabel("Error");
    legend("J(2,3) and Error Pairs", "Nima's J(2,3)", "New J(2,3)")
    title(num2str(trial));
    pause(2);
    clf;
end