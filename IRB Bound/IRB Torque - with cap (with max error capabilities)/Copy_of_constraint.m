%% Constraints function 

%Inputs:

%P - 2x1 vector of impulses (Pn and Pt)
%M - generalized mass matrix
%J - jacobian
%v_pre - pre impact velocity
%v_post - post impact velocity 

%Outputs
%c - constraint equations, written in <= 0 form


function [c, ceq] = constraint(P, M, J, maxWidth, v_pre, v_post)
    %impact law prediction (ILP) for post impact velocity
    ILP = v_pre + inv(M) * J' * [P(1:2)'; P(3) * P(2)]; 
    fake  = v_pre + inv(M) * J(1:2, :)' * P(1:2)'; 
    %initial energy from pre impact velocity 
    initialEnergy = 1/2 * v_pre' * M * v_pre;
    disp("Initial Energy: "+ initialEnergy);
    disp("Final Energy: " + 1/2 * ILP' * M * ILP);
    disp("Final Energy If You Didn't Account for The Additional Torque: " + 1/2 * fake' * M * fake)
    
    %Simplified Energy Ellipse Equaion, written in terms of <= 0
    c(1) = 1/2 * ILP' * M * ILP - initialEnergy; 
    %max width constraint, written in terms of <=
    c(2) = P(3) - maxWidth;
    c(3) = -maxWidth - P(3);
    ceq = [];

end
