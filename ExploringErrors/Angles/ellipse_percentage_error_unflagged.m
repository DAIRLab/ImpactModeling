% Ellipse percentage error unflagged

% Note: in order to run this code, a quick modification must be made in the
% associated function, ellipse_visual. Line 81, which begins with p =
% line... must be commented, as well as lines 84-88 included
%load data

load('ellipse_uniform.mat');

%trial number to draw
desiredTrials = 2000;
for h = 1:desiredTrials
    if sum(bounce_array(h).flags) == 0
        vec(h) = 1;
    else
        vec(h) = 0;
    end
end     
nonflagged = find(vec == 1);
numTrials = sum(vec == 1);
xvec_per = zeros(1,numTrials);
yvec_per = zeros(1,numTrials);

for tr = 1:numTrials
    
trial = nonflagged(tr);

%set up pre-impact ellipse
ang = bounce_array(trial).states(3);
x0 = bounce_array(trial).states(1);
y0 = bounce_array(trial).states(2);
d = (bounce_array(trial).d);   %tangential
n = (bounce_array(trial).n); 
J = [n;d];
x1 = J(1,3)+x0;
y1 = y0 - J(2,3);

[xvec,yvec] = ellipse_visual(ang,x0,y0);
pos = find(yvec == min(yvec));
yper = (y1-min(yvec))/min(yvec)*100;
xper = (x1 - xvec(pos))/xvec(pos)*100;

xvec_per(trial) = xper;
yvec_per(trial) = yper;
disp(tr)
end

M = zeros(2,2);
M(1,1) = mean(xvec_per);
M(1,2) = mean(yvec_per);
M(2,1) = std(xvec_per);
M(2,2) = std(yvec_per);
