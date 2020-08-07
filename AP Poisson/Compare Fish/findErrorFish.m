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


function out = findErrorFish(mu, epsilon, n, d, M, v_pre, v_post)
    predicted = APPoisson_juniors(M, n, d, v_pre, mu, epsilon);
    %Using Impact Law Find Predicted Velocity
    %Find normalized Error Given observed and predicted velocities
    out = norm(v_post - predicted')/norm(v_post);
end