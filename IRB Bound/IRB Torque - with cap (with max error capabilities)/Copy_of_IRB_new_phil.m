%% Newer IRB
clear; 

load('squareDataPhil.mat');

sl = 0.06; %side length of square from data README
m1 = 0.0485; %WHAT IS THE MASS????? Right now we have an area ratio 
                                    %approximation from ellipse

I1 = m1 * sl^2 / 6; %Moment of Inertia

Mass = [m1, 0, 0;
        0, m1, 0;
        0, 0, I1]; %generlaized Mass matrix
    

widths = 0.02;%linspace(0, 0.02, 20);
numTrials = length(squareDataPhil); %Number of Trials
ran = randi([1 2000], 1, numTrials);

count = 0; %counter variable
bestWidth = [];
for a = 1:length(widths) 
    %choose the current max allowable width
    maxWidth = widths(a);
    %reset the error vector matrix
    errorVector = [];
    count = 0;
    for i = 1:length(ran)
        trial = i; %ran(i);
        if true %sum(squareDataPhil(trial).flags) < 1
            % Finding pre and post impact velocities / states
            pre = squareDataPhil(trial).states(4:6)';
            post = squareDataPhil(trial).states(10:12)';

            d = (squareDataPhil(trial).d);   %tangential
            n = (squareDataPhil(trial).n);   %normal

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
            %iterator for non flagged trials 
            count = count + 1;
            %save some data for correlation plots
            error = findError(P, Mass, J, pre, post);
            errorVec(a, count) = error;

            predicted = pre + inv(Mass) * J' * [P(1:2)'; P(3) * P(2)];
            
            bestWidth(1, count) = P(1);
            bestWidth(2,count) = P(2);
            bestWidth(3,count) = P(3);
            useful(1, count) = wrapTo180(rad2deg(squareDataPhil(trial).states(3)));

         end

    end
    
    %find the average error across all trials for a certain max width
end
avErr = mean(errorVec, 2);

figure()
plot(useful, bestWidth(3,:), '.');
xlabel("Wrapped Pre-Impact Angle");
ylabel("Optimal Width [m]")
%% 
figure()
plot(widths, avErr)
%title("IRB with Max Width, Square Data")
xlabel("Max Width [m]");
ylabel("Average Error Across all Trials")