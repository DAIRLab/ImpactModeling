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

function out = findErrorSquare(P, M, J, v_pre, v_post)
    %Using Impact Law Find Predicted Velocity
    predicted = v_pre + inv(M) * J' * P';
    
    s = 0.06; %side length of square
    I = M(3, 3); % Moment of Inertia
    A = s^2; %Area 
    Rg = sqrt(I / A); %radius of gyration
    
    
    %Find normalized Error Given observed and predicted velocities with 
    %omega scaled by the radius of gyration
    v_post_scaled = [v_post(1:2); Rg * v_post(3)];
    predicted_scaled = [predicted(1:2); Rg * predicted(3)];
    out = norm(v_post_scaled - predicted_scaled)/norm(v_post_scaled);
end