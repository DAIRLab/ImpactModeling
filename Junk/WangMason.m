%% Wang Mason Impact Model
%Authors: Andy, Joah, Natalie, Phil

%Goal: 
%Implement the Wang Mason Impact Model from their Paper "Two-Dimensional 
%Rigid-Body Collisions With Friction". This function takes in the 
%pre-impact state of an object and returns the post-impact velocity.

%Inputs:

    %FIRST OBJECT, column vector:
    %m1             mass 
    %p1             radius of gyration of interia
    %x1             x coordinate of center of mass
    %y1             y coordinate of center of mass
    %xdot_1o        pre-impact x velocity 
    %ydot_1o        pre-impact y velocity
    %thetadot_1o    initial rotational velocity 
    
    %SECOND OBJECT, same format as above
    
%Outputs:
    %xdot_1           Post-impact x velocity
    %ydot_1          post-impact y velocity
    %-------
    
%Additional variables used:
    %Px             tangential impulse
    %Py             normal impulse
    %e              coefficient of restitution
    %mu             frictional coefficient

function [output] = WangMason(firstObj, secondObj)
    %% Assign variable names for input vectors
    
    m1 = firstObj(1);
    p1 = firstObj(2);
    x1 = firstObj(3);
    y1 = firstObj(4);
    xdot_1o = firstObj(5);
    ydot_1o = firstObj(6);
    thetadot_1o = firstObj(7);
    
    m2 = secondObj(1);
    p2 = secondObj(2);
    x2 = secondObj(3);
    y2 = secondObj(4);
    xdot_2o = secondObj(5);
    ydot_2o = secondObj(6);
    thetadot_2o = secondObj(7);
    
    
    %% Kinematics: 
    
    %First Object
    %apply angular and impule and momentum laws
    
    %find velocity at contanct points
    xdot_1c = xdot_1 + thetadot_1 * y1; %(8)
    ydot_1c = ydot_1 - thetadot_1 * x1; %(9)

    %Second Object
    %apply angular and impule and momentum laws
    
    %find velocity at contanct points
    xdot_2c = xdot_2 + thetadot_2 * y2; %(13)
    ydot_2c = ydot_2 - thetadot_2 * x2; %(14)
    
    %Determine sliding velocity 
    S = xdot_1c - xdot_2c; %(15)
    
    %Determine compression velocity
    C = ydot_1c - ydot_2c; %(16)
    
    %Solve for constants needed in kinematic equations
    B1 = 1/m1 + 1/m2 + y1^2 / (m1 * p1^2) + y2^2 / (m2 * p2^2); %(19)
    B2 = 1/m1 + 1/m2 + x1^2 / (m1 * p1^2) + x2^2 / (m2 * p2^2); %(20)
    B3 = x1 * y1 / (m1 * p1^2) + x2 * y2 / (m2 * p2^2);         %(21)
    
    S_o = (xdot_1o + thetadot_1o * y1) ...
        - (xdot_2o + thetadot_2o * y2); %(22) 
    C_o = (ydot_1o + thetadot_1o * x1) ...
        - (ydot_2o + thetadot_2o * x2); %(23) 
    
   %% Determine mu
   
   
   
   
   %% Determine e to minimize metric error
   
   
   
   %% Determine Contact Mode 
    
    
    
    
    
    %% Solve for Output
    xdot_1 = (Px/m1)
    

end
