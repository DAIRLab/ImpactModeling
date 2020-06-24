% run ellipse for pre and modified angle

%load data
load('ellipse_uniform.mat');

%trial number to draw
trial = 802;

%specify the angle of the second ellipse
ang2 = -pi/3;

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
[xvec,yvec] = ellipse_visual(ang,x0,y0,'b');
hold on
pos = find(yvec == min(yvec));
disp(min(yvec));
disp(xvec(pos));
plot(xvec(pos),min(yvec),'*','MarkerSize',10,'LineWidth',2,'Color','g')
yper = (y1-min(yvec))/min(yvec)*100;
xper = (x1 - xvec(pos))/xvec(pos)*100;

%draw the second ellipse
hold on
[xvec2,yvec2] = ellipse_visual(ang2,x0,y0,'c');
hold on
pos2 = find(yvec2 == min(yvec2));
plot(xvec2(pos2),min(yvec2),'*','MarkerSize',10,'LineWidth',2,'Color','m')

%plotting the 'ground'
hold on
plot([xlim],[min(yvec2) min(yvec2)],'LineWidth',2,'Color','k')
hold on
plot([xlim],[y1 y1],'LineWidth',2,'Color','k')

title("Trial #: " + trial, 'FontSize',15);
xlabel('X-axis','FontSize',15);
ylabel('Y-axis','FontSize',15);
%annotation('textbox',[.13 0.7 .6 .2],'String',['Jacobian contact x and y = ' num2str(x1) ' and ' num2str(y1)],'EdgeColor','none','FontSize',10)
%annotation('textbox',[.13 0.665 .6 .2],'String',['Ellipse contact x and y = ' num2str(xvec(pos)) ' and ' num2str(min(yvec))],'EdgeColor','none','FontSize',10)
%annotation('textbox',[.13 0.63 .6 .2],'String',['Percentage change = ' num2str(xper) '% and ' num2str(yper) '%'],'EdgeColor','none','FontSize',10)
legend('Jacobian contact','Ellipse','Ellipse contact','Ellipse2','Ellipse2 contact')
