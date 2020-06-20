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
    

widths = linspace(0, 0.01, 20);
numTrials = 100; %Number of Trials
ran = randi([1 2000], 1, numTrials);

count = 0; %counter variable

for a = 1:length(widths) 
    %choose the current max allowable width
    maxWidth = widths(a);
    %reset the error vector matrix
    errorVector = [];
    
    for i = 1:length(ran)
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
            
            %add the trial's error to error vector
            count = count + 1;
            errorVector(count) = findError(P, Mass, J, pre, post);

         end

    end
    
    %find the average error across all trials for a certain max width
    avgError(a) = mean(errorVector);
end

% disp("Tangential Impulse: " + P(1) + " [N*s]")
% disp("Normal Impluse: " + P(2) + " [N*s]")

%% Post Process Data
plot(widths, avgError)
%title("Maximum Allowable Width vs. Error Plot")
xlabel("Maximum Allowable Width [m]")
ylabel("Average Normalized Error Across All Trials")

