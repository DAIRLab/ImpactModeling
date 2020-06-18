%% run ellipse for pre and post data to visualize impacts

%load data
load('ellipse_uniform.mat');

%trial number to draw
trial = 13;

%set up pre-impact ellipse
ra = 0.035;
rb = 0.025;
ang = bounce_array(trial).states(3);
x0 = bounce_array(trial).states(1);
y0 = bounce_array(trial).states(2);

%create a new figure window 
figure()

%draw pre-impact
ellipse(ra, rb, ang, x0, y0, 'b')

%set up post-impact ellipse
ang1 = bounce_array(trial).states(9);
x1 = bounce_array(trial).states(7);
y1 = bounce_array(trial).states(8);

%draw post-impact 
ellipse(ra, rb, ang1, x1, y1, 'r')
hold on

title("Trial #: " + trial);
hold off