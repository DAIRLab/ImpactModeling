syms x y h k a b th

eq1 = (((x-h)*cos(th) + (y-k)*sin(th))^2)/a^2 + (((x-h)*sin(th) - (y-k)*cos(th))^2)/b^2 == 1;
eq1_2 = (((x)*cos(th) + (y)*sin(th))^2)/a^2 + (((x)*sin(th) - (y)*cos(th))^2)/b^2 == 1;

x1 = solve(eq1, sym('x'));
x1_2 = solve(eq1_2, sym('x'));