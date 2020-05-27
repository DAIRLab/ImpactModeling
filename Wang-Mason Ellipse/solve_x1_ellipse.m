  
function x1 = solve_x1_ellipse(y1, h0, k0, theta1) 

syms x y h k a b th

eq1 = (((x-h)*cos(th) + (y-k)*sin(th))^2)/a^2 + (((x-h)*sin(th) - (y-k)*cos(th))^2)/b^2 == 1;
x1eq = solve(eq1, sym('x'));


a0 = 0.07/2; %semi axis major
b0 = 0.05/2; %semi axis minor
h0 = 0; %h0 and k0 are 0 since we assume the com is the origin of the ellipse
k0 = 0;
th0 = theta1; %tilting angle
y0 = y1; %distance from com to contact point (going down from com to contact point hence the negative sign) because y1 
%goes up from the table ie contact point to COM (

%solving for x1v
x1 = subs(x1eq, [a,b,h,k,th,y], [a0,b0,h0,k0,th0,y0]); %gives out the x distance from the COM to the contact point in both positive
%negative since the answer is a square root
