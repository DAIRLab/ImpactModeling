%% Wang Mason Model

%Takes in the pre-impact state and returns the predicted post-impact
%state using the Wang Mason Model

%Inputs
%   m - objects mass
%   S_0 - Initial sliding velocity
%   C_0 - Initial compression velocity
%   B - 3x1 column vector containing B1, B2, and B3 (constants)
%   v_0 - 3x1 column vector with initial x, y, and rotational velocity
%   u - mu, coefficiant of friction 
%   e - epsilon, coefficiant of restitution


%Outputs
%   v_1 - 2x1 column vector with predicted post impact x and y velocities

function [v_1]=wang_juniors(m, S_0, C_0, B, v_0, u ,e)
    %Assign B vector to propper B variables
    B1 = B(1);
    B2 = B(2);
    B3 = B(3);
    
    u_s = -B3/B1; %(34)
    %Initial Sign of Sliding Velocity 
    if (S_0 ~= 0) 
        s = sign(S_0);
    else
        s = 1;
    end


    %Use Table 1 to determine modes (conditionals)
    Pd = (B2 + s * u * B3) * s * S_0;  %(35)
    Pq = (u * B1 + s * B3)*(-C_0);     %(36)
    %Apply equations 39 - 48 based on mode
    %Sliding (Second Row of Table)
    
    if (Pd > (1+e)*Pq)
        Py = - (1+e) * C_0 / (B2 + s * u * B3); %(40)
        Px = - s * u * Py;                      %(39)

        %R (Third Row of Table)
    elseif (Pq < Pd) && (Pd < (1+e) * Pq)
        disp("R")
        if u > abs(u_s) %R-Sticking
            Py = -(1+e)* C_0 /(B2 + s * u * B3);  %(44)
            Px = (B3*Py - S_0) / B1;              %(43)

        else %R-Reversed Sliding
            Py = -(1+e) * C_0/(B2 + s * u * B3);              %(48)
            Px = s * u * (Py - 2 * S_0 / (B3 + s * u * B1));  %(47)
        end

    %C (Fourth Row of Table)
    elseif (Pd < Pq)
        if u > abs(u_s) %C-Sticking
             disp("C-stick");
             Py = -(1+e) * (B1 * C_0 + B3 * S_0)/(B1*B2 - B3^2); %(42)
             Px = (B3 * Py - S_0)/B1;                            %(41)
        else %C-Reversed Sliding
             disp("C-slide");
             Py = -(1+e)/(B2 - s * u * B3)* (C_0 + ...
                 (2*s*u*B3*S_0)/(B3 + s * u * B1));               %(46)
             Px = s * u *(Py - 2 * S_0 / (B3 + s * u * B1));    %(45)
         end    
    else
        disp("Error, none of contact mode requirements met");
    end    

    % calculate the post impact velocities according to the contact mode
    x1dot_calc = Px/m + v_0(1); 
    y1dot_calc = Py/m + v_0(2); 
    
    %format into output vector
    v_1 = [x1dot_calc; y1dot_calc];
        
        
end
    