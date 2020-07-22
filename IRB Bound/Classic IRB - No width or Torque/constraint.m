%% Constraints function 

%Inputs:

%P - 2x1 vector of impulses (Pn and Pt)
%M - generalized mass matrix
%J - jacobian
%v_pre - pre impact velocity
%v_post - post impact velocity 

%Outputs
%c - constraint equations, written in <= 0 form


function [c, ceq] = constraint(P, M, J, v_pre, v_post)
    %impact law prediction (ILP) for post impact velocity
    ILP = v_pre + inv(M) * J' * P'; 
    %initial energy from pre impact velocity 
    initialEnergy = 1/2 * v_pre' * M * v_pre;
    %Simplified Energy Ellipse Equaion, written in terms of <= 0
    c(1) = 1/2 * ILP' * M * ILP - initialEnergy; 
    
    ceq = [];

end
