load('ellipse_uniform.mat');

n = 190;
pre = bounce_array(n).states(1:6); 
post = bounce_array(n).states(7:12);
J = [bounce_array(n).d; bounce_array(n).n];
%J = [1 0 rx; 0 1 ry]
u = 0.8;
e = 0.8;

%[stick,Mu,Ep] = ErrorEllipse(n,pre,post,J)
%[stick, errors,bestMu, bestEpsilon] = ErrorEllipseFriday(n,pre,post,J);
[error] = ErrorEllipseFridaySingle(pre,post,J,u,e)

%x = 0.01:0.01:.99;
%y = 0.01:0.01:.99;
%contourf(x,y,errors,20)
%colorbar
