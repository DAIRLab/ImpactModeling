
function [stick, bestMu, bestEpsilon] = ErrorWang(n,pre,post,J)

%Establishing Variables
stepSize = 0.01;
a0 = 0.07/2; 
b0 = 0.05/2; 
rho = 0.5 * sqrt(a0^2 + b0^2); 
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4; 

%Initialize Error Matrix
sz = 1/stepSize - 1;
errors = zeros(sz,sz); 

%Pre Impact Conditions (object 1 COM)
x1_0 = pre(1,1);
y1_0 = pre(1,2); 
th1_0 = pre(1,3);
x1dot_0 = pre(1,4);
y1dot_0 = pre(1,5);
th1dot_0 = pre(1,6);

%Pre Impact Conditions (object 2 COM - arbitrary)
x2_0 = 1; 
y2_0 = 1;
th_0 = 1;
x2dot_0 = 0;
y2dot_0 = 0;
th2dot_0 = 0;

% signx1 = 0;
% if th1_0 > 0 %turning anticlockwise
%     r = rem(th1_0,pi)/pi;
%     if r > 0 && r < 0.25 || r > 0.5 && r < 0.75 %angle between 0 to 90 degrees and 180-270
%         signx1 = 1; %x1 is to the right of contact point
%     elseif r > 0.25 && r < 0.5 || r > 0.75 %angle between 90 and 180 and 270-360
%         signx1 = -1; %x1 will be on the left
%     else %when it is right on top
%         signx1 = 1;
%     end
% elseif th1_0 < 0 %turning clockwise
%     r = rem(abs(th1_0),pi)/pi;
%     if r > 0 && r < 0.25 || r > 0.5 && r < 0.75 %angle between 0 to 90 degrees and 180-270
%         signx1 = -1; %x1 is to the left of contact point
%     elseif r > 0.25 && r < 0.5 || r > 0.75 %angle between 90 and 180 and 270-360
%         signx1 = 1; %x1 will be on the right
%     else %when it is right on top
%         signx1 = 1;
%     end
% elseif th1_0 == 0 %doesn't even rotate
%     signx1 = 1;
% else
%     signx1 = 1; 
% end


%Finding the contact points
a = a0;
b = b0;
th = th1_0;
k = 0;
h = 0;
y = -y1_0;

%x1d = sqrt(0.06^2/2 - pre(1,2)^2);
x1c = double((a*b*(b^2*cos(th)^2 - k^2*cos(th)^4 + a^2*sin(th)^2 - y^2*cos(th)^4 - k^2*sin(th)^4 - y^2*sin(th)^4 + 2*k*y*cos(th)^4 + 2*k*y*sin(th)^4 - 2*k^2*cos(th)^2*sin(th)^2 - 2*y^2*cos(th)^2*sin(th)^2 + 4*k*y*cos(th)^2*sin(th)^2)^(1/2) + b^2*h*cos(th)^2 + a^2*h*sin(th)^2 - a^2*k*cos(th)*sin(th) + b^2*k*cos(th)*sin(th) + a^2*y*cos(th)*sin(th) - b^2*y*cos(th)*sin(th))/(a^2*sin(th)^2 + b^2*cos(th)^2));

%x1c = signx1*abs(x1c);

%Soving for constants using Jacobian
V_c = J*[pre(1,4);pre(1,5);pre(1,6)];
S_0 = V_c(2);
C_0 = V_c(1);

%Finding B1, B2, and B3
B1 = 1 + y/(m1*rho^2); 
B2 = 1 + x1c/(m1*rho^2); 
B3 = x1c * y/m1/rho^2;    
%Static friction coefficiant 
u_s = -B3/B1; %(34)
%Initial Sign of Sliding Velocity 
if (S_0 ~= 0) 
    s = sign(S_0);
else
    s = 1;
end

%Double for loop
for a = 1:sz %varying mu from [0, 1]
    u = stepSize * a;
    
    Pd = (B2 + s * u * B3) * s * S_0;  %(35)
    Pq = (u * B1 + s * B3)*(-C_0);     %(36)
    
    for b = 1:sz 
        %Use Table 1 to determine modes (conditionals)
        %Apply equations 39 - 48 based on mode
        e = (stepSize) * b ;
        
        %Second Row
        % 1. Sliding
        if (Pd > (1+e)*Pq)
            Py = - (1+e) * C_0 / (B2 + s * u * B3); %(40)
            Px = - s * u * Py;                      %(39)
        
        %Third Row
        elseif (Pq < Pd) && (Pd < (1+e)*Pq)
            % 2. R-sticking
            if u > abs(u_s) 
                Py = -(1+e)* C_0 /(B2 + s * u * B3);  %(44)
                Px = (B3*Py - S_0) / B1;              %(43)
            % 3. R-Reversed Sliding
            else 
                Py = -(1+e) * C_0/(B2 + s * u * B3);              %(48)
                Px = s * u * (Py - 2 * S_0 / (B3 + s * u * B1));  %(47)
            end
            
        %Fourth Row
        elseif (Pd < Pq)
            % 4. C-Sticking
            if u > abs(u_s) 
                 Py = -(1+e) * (B1 * C_0 + B3 * S_0)/(B1*B2 - B3^2); %(42)
                 Px = (B3 * Py - S_0)/B1;                            %(41)
            % 5. C-Reversed Sliding
            else 
                 Py = -((1+e)/(B2 - s * u * B3)) * (C_0 + ...
                     (2*s*B3*S_0)/(B3 + s * u * B1));               %(46)
                 Px = s * u *(Py - (2 * S_0/(B3 + s * u * B1)));    %(45)
             end    
        else
            disp("Error, none of contact mode requirements met");
        end    

        % calculate the post impact velocities according to the contact
        % mode
        x1dot_calc = Px/m1 + pre(1,4); 
        y1dot_calc = Py/m1 + pre(1,5); 
        V_cpost = J*[post(1,4);post(1,5);post(1,6)];
        x_post = V_cpost(2);
        y_post = V_cpost(1);

        % calculate error via least squares m
        error = sqrt((x_post - x1dot_calc)^2 + (y_post - y1dot_calc)^2)/sqrt(y_post^2 + x_post^2); 
        % input error into error matrix
        errors(a, b) = error;
    end
end

% determine minimum error 
minimum = min(min(errors));
% determine the indices of the minimum error
[i,j] = find(errors == minimum);
% determine the "best" mu and epsilon value which yields the minimum error
bMu = i * stepSize;
bEpsilon = j * stepSize;
    
    [g,h] = size(bMu);
    if g>1 || h>1 
        stick = 1; %multiple best mus - sticking
    else
        stick = 0; %a single best my - no sticking
    end
    
    bestMu = min(bMu);
    bestEpsilon = min(bEpsilon);