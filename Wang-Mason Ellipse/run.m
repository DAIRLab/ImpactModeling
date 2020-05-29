load('ellipse_uniform.mat');

n = 40;
pre = bounce_array(n).states(1:6); 
post = bounce_array(n).states(7:12);
J = [bounce_array(n).n; bounce_array(n).d];

%[stick,Mu,Ep] = ErrorEllipse(n,pre,post,J)
[stick, bestMu, bestEpsilon] = ErrorEllipse(n,pre,post,J)