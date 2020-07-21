function [P,w,error] = IRB_NEW_torque_square(Pre,Post,trial,p1,p2)

% Set up Constants
sl = 0.06; %side length of square from data README
rho = sqrt(sl^2/6); %using I/m where I = m *s^4 / 12
m1 = 1; %cancels out as explained in variables section above
I1 = m1 * sl^2 / 6; % moment of inertia of square
Mass = [m1, 0, 0;
        0, m1, 0; 
        0, 0, I1]; 
    
d = [1 0 Pre(2)];
n = [0 1 Pre(1)];

n(3) = n(3)*p1;
d(3) = d(3)*p2;

J = [d;n;0 0 1]; %Jacobian

pre = Pre(4:6)';
post = Post(4:6)';

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
