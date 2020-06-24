%{
This script will run the IRB_NEW_torque function a specified number of
times and compute the average width of the contact point.
%}
load('ellipse_uniform.mat')

numTrials = 1000;
ran = randi([1, 2000], 1, numTrials);
avgW = 0;
avgdx = 0;
avgdy = 0;

for i=1:numTrials
    z = ran(i);
    [vector, w] = IRB_NEW_torque(z);
    avgW = avgW + w;
    
    avgdx = avgdx + 10*bounce_array(z).states(1);
    avgdy = avgdy + 20*bounce_array(z).states(2);
end

avgW = avgW / numTrials
avgdx = avgdx / numTrials
avgdy = avgdy / numTrials
