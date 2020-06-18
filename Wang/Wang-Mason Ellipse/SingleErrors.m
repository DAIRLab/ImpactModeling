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
    
figure

pd = fitdist(errors','Normal');
step = 2/(length(errors)-1);
x = 0:step:2;
y = pdf(pd,x);
plot(x,y,'LineWidth',2)
standardEp = pd.sigma;
yl = ylim;
title({'Density Function for 2000 Cases', 'Mean Error = ' num2str(mean(errors))})
xlabel('Normalized Error Metric')
ylabel('Probability Density Estimate')
    
