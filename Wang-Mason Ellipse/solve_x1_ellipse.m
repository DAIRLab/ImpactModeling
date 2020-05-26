function x1 = solve_x1_ellipse

syms x y h k a b th

eq1 = (((x-h)*cos(th) + (y-k)*sin(th))^2)/a^2 + (((x-h)*sin(th) - (y-k)*cos(th))^2)/b^2 == 1;
x1eq = solve(eq1, sym('x'));


a0 = 0.07/2; %semi axis major
b0 = 0.05/2; %semi axis minor
h0 = 0; %h0 and k0 are 0 since we assume the com is the origin
k0 = 0;
th0 = theta1; %tilting angle
y0= y1; %distance from com to contact point

%solving for x1v
x1 = subs(x1eq, [a,b,h,k,th,y], [a0,b0,h0,k0,th0,y0]); %gives out the x distance from the COM to the contact point in both positive
%negative since the answer is a square root

 
