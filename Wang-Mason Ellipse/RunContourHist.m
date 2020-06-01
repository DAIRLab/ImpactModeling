load('ellipse_uniform.mat');
errorssum = zeros(99,99);
muvec = 0.01:0.01:(1-0.01);
epvec = 0.01:0.01:(1-0.01);
minerrvec = zeros(1,2000);

for n = 1:2000
    
    pre = bounce_array(n).states(1:6);
    post = bounce_array(n).states(7:12);
    J = [bounce_array(n).n; bounce_array(n).d];
    
    %[stick,Mu,Ep] = ErrorEllipse(n,pre,post,J)
    [stick, bestMu, bestEpsilon,errors,minerr] = ErrorEllipseFriday(n,pre,post,J);
    errorssum= errorssum +errors;
    minerrvec(n)= minerr;
    
end

% errorsavg = errorssum/2000;
% 
% contourf(epvec,muvec,errorsavg,40)
histogram(minerrvec,'Normalization','Probability');