function [P,w,error] = IRB_NEW_torque(bounce_array,trial,p1,p2)

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
Mass = [m1, 0, 0;
        0, m1, 0; 
        0, 0, I1]; 

% Finding pre and post impact velocities / states
pre = bounce_array(trial).states(4:6)';
post = bounce_array(trial).states(10:12)';

d = (bounce_array(trial).d);   %tangential
n = (bounce_array(trial).n);   %normal
n(3) = n(3)*p1;
d(3) = d(3)*p2;

J = [d;n;0 0 1]; %Jacobian

fun = @(P)(findError_torque3(P, Mass, J, pre, post));
nonlcon = @(P)(constraint_torque(P, Mass, J, pre, post));

P0 = [0 0 0];
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
options = optimoptions('fmincon','FiniteDifferenceType','central', ...
                    'StepTolerance',1e-15, 'Display','off');
P = fmincon(fun, P0, A, b, Aeq, beq, lb, ub, nonlcon);

%w = abs(P(3))/P(2);
w = P(3)/P(2);

J2 = [d;n];
P2 = P(1:2);

predicted = pre + inv(Mass) * J' * P';
%error = (norm(post - predicted)/norm(post))^2;
%error = norm(post - predicted);
error = (post(1:2) - predicted(1:2))'*(post(1:2) - predicted(1:2))

% disp("Tangential Impulse: " + P(1) + " [N*s]")
% disp("Normal Impulse: " + P(2) + " [N*s]")
% disp("Torque: " + P(3) + "[N*m]")
% disp("w: " + w + "[m]")

end