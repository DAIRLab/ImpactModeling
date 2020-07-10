 
% Find The Error Given P

%Inputs
%P - Impulse vector (3x1 with Pt, Pn, and Torque)
%M - generalized mass matrix
%J - jacobian
%v_pre - pre impact velocity
%v_post - observed post impact velocity

%Outputs
%out = Normalized error of observed post impact velocity and predicted post
%impact velocity using the impact law with the P vector

function out = findError_torque3(P, M, J, v_pre, v_post)
    %Using Impact Law Find Predicted Velocity
    predicted = v_pre + inv(M) * J' * P';
    %Find normalized Error Given observed and predicted velocities 
    out = (v_post(1:2) - predicted(1:2))'*(v_post(1:2) - predicted(1:2));
end