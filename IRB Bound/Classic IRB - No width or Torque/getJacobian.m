%% Get Jacobian for Ellipse
%Inputs:
%Angle, COM Position
%Outputs:
%Contact Jacobian 

function J = getJacobian(vec)
    x = vec(1);
    y = vec(2);
    ang = vec(3);

    [xvec,yvec] = ellipse_visual(vec); %find all points on ellipse
    
    idx = find(min(yvec)==yvec);  %find index of mimimum y point (contact point)
    contactPoint = [xvec(idx), yvec(idx)]; %find contact point
    %Construct Jacobian
    n = [0, 1, contactPoint(1)-x];
    d = [1, 0, y - contactPoint(2)];
    J = [d ; n];

end