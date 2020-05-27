% This function calculates the coefficient of restitution directly from the
% ellipse data given, since the data provides pre and post impact
% velocities.

function avgEps = calcRestitution

load('ellipse_uniform.mat');

% initialize vector to hold calculated coefficient of restitution values
calcEps = zeros(200,1);
ran = randi([1 2000],1,200);
% iterate through data 
for i = 1:200
    n = ran(i);
    %accounting for x and y components of velocity
    preVel = sqrt(bounce_array(n).states(4)^2 + bounce_array(n).states(5)^2);
    postVel = sqrt(bounce_array(n).states(10)^2 + bounce_array(n).states(11)^2);
   
    coef = postVel / preVel;
    
    calcEps(i) = coef;
end

avgEps = sum(calcEps) / 200;

end


