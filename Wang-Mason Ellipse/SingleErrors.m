load('ellipse_uniform.mat');

e = 0.6385;
u = 0.1624;

errors = zeros(1,2000);

for n = 1:2000
    pre = bounce_array(n).states(1:6); 
    post = bounce_array(n).states(7:12);
    J = [bounce_array(n).n; bounce_array(n).d];
    [error] = ErrorEllipseFridaySingle(pre,post,J,u,e);
    errors(n) = error;
end

histogram(errors)
    