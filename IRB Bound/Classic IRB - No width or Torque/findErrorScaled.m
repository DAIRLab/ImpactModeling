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

function out = findErrorScaled(P, M, J, v_pre, v_post, scale)
    %Using Impact Law Find Predicted Velocity
    predicted = v_pre + inv(M) * J' * P';
    %Find normalized Error Given observed and predicted velocities 
    out = (norm(v_post(1:2) - predicted(1:2))/norm(v_post(1:2)))^2 + ... 
    (norm(v_post(3) - predicted(3))/norm(v_post(3)))^2 / scale;
end