function x1v = solve_x1_ellipse(y0,h0,k0,th0)

syms x y h k a b th

eq1 = (((x-h)*cos(th) + (y-k)*sin(th))^2)/a^2 + (((x-h)*sin(th) - (y-k)*cos(th))^2)/b^2 == 1;
eq1_2 = (((x)*cos(th) + (y)*sin(th))^2)/a^2 + (((x)*sin(th) - (y)*cos(th))^2)/b^2 == 1;

x1 = solve(eq1, sym('x'));
x1_2 = solve(eq1_2, sym('x'));

%solving for x1v
x1v = subs(x1, [a,b,h,k,th,y], [0.07,0.05,h0,k0,th0,y0]);

