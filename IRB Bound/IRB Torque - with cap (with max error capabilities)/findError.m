%% Find The Error Given P

%Inputs
%P - Impulse vector (2x1 with Pn and Pt)
%M - generalized mass matrix
%J - jacobian
%v_pre - pre impact velocity
%v_post - observed post impact velocity

%Outputs
%out = Normalized error of observed post impact velocity and predicted post
%impact velocity using the impact law with the P vector

function out = findError(P, M, J, v_pre, v_post)
    %Using Impact Law Find Predicted Velocity
    predicted = v_pre + inv(M) * J' * [P(1:2)'; P(3) * P(2)];
    %Find normalized Error Given observed and predicted velocities 
    out = norm(v_post - predicted)/norm(v_post);
end
