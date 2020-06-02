load('ellipse_uniform.mat');
errorssum = zeros(99,99);
muvec = 0.01:0.01:(1-0.01);
epvec = 0.01:0.01:(1-0.01);
minerrvec = zeros(1,2000);
count =0;

for n = 1:2000
    
    pre = bounce_array(n).states(1:6);
    post = bounce_array(n).states(7:12);
    J = [bounce_array(n).n; bounce_array(n).d];
    
    [stick, bestMu, bestEpsilon,errors,minerr] = ErrorEllipseFriday(n,pre,post,J);
    errorssum= errorssum +errors;
    minerrvec(n)= minerr;
    
    %contour plot for every trial
    %     contourf(epvec,muvec,errors,40)
    %     pause(2)
    
    
    %contour plot for every sliding trial
    %     if stick ==1
    %         count = count+1;
    %     else
    %
    %         contourf(epvec,muvec,errors,200)
    %         ylim([0 0.5])
    %         pause(2)
    %     end
    %
end

errorsavg = errorssum/2000;

contourf(epvec,muvec,errorsavg,250)
colorbar
%histogram(minerrvec,'Normalization','Probability');
