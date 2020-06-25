% drawing ellipses

%load data
load('ellipse_uniform.mat');

%trial number to draw
trial = 200;

%set up pre-impact ellipse
ang = bounce_array(trial).states(3);
x0 = bounce_array(trial).states(1);
y0 = bounce_array(trial).states(2);
d = (bounce_array(trial).d);   %tangential
n = (bounce_array(trial).n);   %normal

%create a new figure window 
figure()

%draw pre-impact 
[xvec,yvec] = ellipse_visual(ang,x0,y0,'b');
hold on
pos = find(yvec == min(yvec));
disp(min(yvec));
disp(xvec(pos));

%plot the variation points
range = ang*0.4;
sample_ang = linspace(ang-range, ang+range, 51);
x1_vec = zeros(1,51);
y1_vec = zeros(1,51);

for i = 1:51
    ang2 = sample_ang(i);
    [xvec3,yvec3] = ellipse_vector(ang2,x0,y0);
    pos3 = find(yvec3 == min(yvec3));
    x1_vec(i) = min(xvec3(pos3));
    y1_vec(i) = min(yvec3);
end

%split between positive and negative y values
a = find(y1_vec >= 0);
b = find(y1_vec < 0);
hold on
plot(x1_vec(a),y1_vec(a),'o','MarkerEdgeColor','g')
hold on
plot(x1_vec(b),y1_vec(b),'o','MarkerEdgeColor','r')
title("Trial #: " + trial, 'FontSize',15);
xlabel('X-axis','FontSize',15);
ylabel('Y-axis','FontSize',15);
hold on
plot([xlim],[0 0],'LineWidth',2,'Color','k')
hold on
plot(xvec(pos),min(yvec),'*','MarkerSize',10,'LineWidth',2,'Color','c') 
hold on
plot(x0+n(1,3),y0-d(1,3),'*','MarkerSize',10,'LineWidth',2,'Color','y') 
legend('Ellipse','Permissible contacts','Penetration contacts','Ground','Configuration contact','Jacobian contact')
