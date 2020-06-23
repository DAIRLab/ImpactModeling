% run ellipse for pre and post data to visualize impacts

%load data
load('ellipse_uniform.mat');

%trial number to draw
trial = 1000;

%set up pre-impact ellipse
ang = bounce_array(trial).states(3);
x0 = bounce_array(trial).states(1);
y0 = bounce_array(trial).states(2);
d = (bounce_array(trial).d);   %tangential
n = (bounce_array(trial).n); 
J = [n;d];
x1 = J(1,3)+x0;
y1 = y0 - J(2,3);

%create a new figure window 
figure()

%draw pre-impact 
plot(x1,y1,'*','MarkerSize',10,'LineWidth',2,'Color','r')
hold on
ellipse_visual(ang, x0, y0)
hold on
plot([xlim],[y1 y1],'LineWidth',2,'Color','k')
title("Trial #: " + trial, 'FontSize',15);
xlabel('X-axis','FontSize',15);
ylabel('Y-axis','FontSize',15);
