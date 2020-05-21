%% Find Contact Mode and Apply Equations
% Use Table 1 to determine modes (conditionals)
% Apply equations 39 - 48 based on mode

%Sliding (Second Row of Table)
if (Pd > (1+e)*Pq)
    Py = - (1+e) * C_0 / (B2 + s * u * B3);
    Px = - s * mu * Py;
%R (Third Row of Table)
elseif (Pq < Pd) && (Pd < (1+e) * Pq)
    
    if mu > abs(mu_s) %R-Sticking
        Py = -(1+e)* C_0 /(B2 + s * u * B3);
        Px = (B3*Py - S_0) / B1;
    
    else %R-Reversed Sliding
        Py = -(1+e) * C_0/(B2 + s * u * B3);
        Px = s * u * (Py - 2 * S_0 / (B3 + s * u * B1)); 
    end
%C (Fourth Row of Table)   
elseif (Pd < Pq)
    
     if mu > abs(mu_s) %C-Sticking
         Py = -(1+e) * (B1 * C_0 + B3 * S_0)/(B1*B2 - B3^2);
         Px = (B3 * Py - S_0)/B1;
     else %C-Reversed Sliding
         Py = -(1+e)/(B2 - s * u * B3)* (C_0 + ...
             (2*s*B3*S_0)/(B3 + s * u * B1));
         Px = s * u *(Py - 2 * S_0 / (B3 + s * u * B1)); 
     end    
else
    disp("Error, none of contact mode requirements met");
end
