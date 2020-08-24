%% Get Jacobian for Ellipse
%Inputs:
%Angle, COM Position
%Outputs:
%Contact Jacobian 

function J = getJacobian(x, y, ang)

    [xvec,yvec] = ellipse_visual(ang, x, y); %find all points on ellipse
    
    idx = find(min(yvec)==yvec);  %find index of mimimum y point (contact point)
    contactPoint = [xvec(idx), yvec(idx)]; %find contact point
    %Construct Jacobian
    n = [0, 1, x - contactPoint(1)];
    d = [1, 0, y - contactPoint(2)];
    J = [d ; n];

end