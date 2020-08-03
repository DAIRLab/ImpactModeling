function Analysis(n)

%ReadData
%Junior Team

%This code takes as an input a csv file that contains information about the
%object's position and velocity in the plane. The code will identify the
%regions of the data where an impact has occurred, and will output the
%pre-impact and post-impact set position and velocity of the object

%Be sure to change this depending on which computer you are using this. Use
%the command window to find the directory of where your files are
D = readmatrix(['/Users/andyeske/Desktop/Simulations/dice-data-processed/traj_' num2str(n) '.csv']);

%Column 1: x position 
x = D(:,1);
%Column 2: z position
z = D(:,2);
%Column 3: theta
th = D(:,3);
%Column 4: xDot
xDot = D(:,4);
%Column 5: zDot
zDot = D(:,5);
%Column 6: thetaDot
thDot = D(:,6);

l = length(zDot);

for jj = 1:l
    
     the = th(jj);
     cc1 = -x(1);
     cc2 = z(1);
     ax = 0.03;
     ay = 0.03;
     bx = 0.03;
     by = - 0.03;
     cx = - 0.03;
     cy = - 0.03;
     dx = - 0.03;
     dy = 0.03;

     a = [ax;ay];
     b = [bx;by];
     c = [cx;cy];
     d = [dx;dy];
     
     vx = [a(1),b(1),c(1),d(1),a(1)];
     vy = [a(2),b(2),c(2),d(2),a(2)];

     R = [cos(the),sin(the);-sin(the),cos(the)];

     a2 = R*a;
     b2 = R*b;
     c2 = R*c;
     d2 = R*d;

     vx = [a2(1),b2(1),c2(1),d2(1),a2(1)]+cc1;
     vy = [a2(2),b2(2),c2(2),d2(2),a2(2)]+cc2;
    
end

%To identify the positions at which an impact occurred, we will take a look
%at the changes in zDot over time - in essence, we will try to
%find the parts where zDot changes sign. When an object hits the ground
%from above, it velocity should reverse. Hence, we implement the
%following code:

imp = [];

%This for loop finds a vector imp, which tracks all the times zDot changes
%sign
for k = 2:l
    if sign(zDot(k)) ~= sign (zDot(k-1))
        imp = [k, imp];
    end  
end

l2 = length(imp);
imp = imp - 1;
impR = [];

%This for loop compares the previous imp vector against the z positions, to
%find only the values at which z is at a local minima
for h = 1:l2
    a = imp(h);
    if z(a) < z(a+1) && z(a) < z(a-1) 
        impR = [a,impR];    
    end
end

%output: two matrices that contains all the information for the pre and
%post impacts
impR2 = impR + 1;
miny = [];
minx = [];

for hh = 1:length(impR)
     
     i = impR(hh);
     the = th(i);
     c1 = -x(i);
     c2 = z(i);
     ax = 0.03;
     ay = 0.03;
     bx = 0.03;
     by = - 0.03;
     cx = - 0.03;
     cy = - 0.03;
     dx = - 0.03;
     dy = 0.03;

     a = [ax;ay];
     b = [bx;by];
     c = [cx;cy];
     d = [dx;dy];
     
     vx = [a(1),b(1),c(1),d(1),a(1)];
     vy = [a(2),b(2),c(2),d(2),a(2)];

     R = [cos(the),sin(the);-sin(the),cos(the)];

     a_2 = R*a;
     b_2 = R*b;
     c_2 = R*c;
     d_2 = R*d;

     vx = [a_2(1),b_2(1),c_2(1),d_2(1),a_2(1)]+c1;
     vy = [a_2(2),b_2(2),c_2(2),d_2(2),a_2(2)]+c2;
     
     miny = [miny, min(vy)];
     aa = find(vy == min(vy));
     minx = [minx, vx(min(aa))];
     
end

Pre = [impR', x(impR), z(impR), th(impR), xDot(impR), zDot(impR), thDot(impR), minx'+x(impR),z(impR)-miny'];
Post = [impR2', x(impR2), z(impR2), th(impR2), xDot(impR2), zDot(impR2), thDot(impR2)];

vect = (impR(1)-5):1:(impR(1)+5)
time = (1/250)*vect;

% Impact plot +- 5 points

plot(-x(vect),z(vect),'r*')
title({['First Impact for trial ', num2str(n)], '+- 5 data points'})
xlabel('x position (m)')
ylabel('z position (m)')
set(gca, 'FontSize', 12)

% x position and xDot vs time
figure
plot(time,-x(vect))
hold on
plot(time,xDot(vect))
title('x position and xDot vs time')
xlabel('Time (s)')
ylabel('Magnitude')
legend('x position', 'x velocity')
set(gca, 'FontSize', 12)

% z position and zDot vs time
figure
plot(time,z(vect))
hold on
plot(time,zDot(vect))
title('z position and zDot vs time')
xlabel('Time (s)')
ylabel('Magnitude')
legend('z position', 'z velocity')
set(gca, 'FontSize', 12)

%z position fitted
figure 
tiledlayout(2,1)
nexttile
p = polyfit(time(1:6),z(vect(1:6))',1)
x1 = time(1):0.0001:time(6);
y = p(1)*x1 + p(2);
p2 = polyfit(time(6:11),z(vect(6:11))',1)
x2 = time(6):0.0001:time(11);
y2 = p2(1)*x2 + p2(2);
plot(time(1:5),z(vect(1:5)),'ro','MarkerSize',8,'MarkerFaceColor','r')
hold on
plot(time(7:11),z(vect(7:11)),'bo','MarkerSize',8,'MarkerFaceColor','b')
hold on
plot(time(6),z(vect(6)),'ko','MarkerSize',8,'MarkerFaceColor','k')
hold on
plot(x1,y,'r')
hold on
plot(x2,y2,'b')
hold on
plot([time(6), time(6)], ylim, 'k--')
legend('Pre-Impact points', 'Post-impact points', 'Impact point','Pre-linear fit','Post-linear fit')
title({['Z position analysis for trial ', num2str(n)], '+- 5 data points'})
xlabel('Time (s)')
ylabel('Position (m)')
set(gca, 'FontSize', 12)

nexttile
plot(time(1:5),zDot(vect(1:5)),'ro','MarkerSize',8,'MarkerFaceColor','r')
hold on
plot(time(7:11),zDot(vect(7:11)),'bo','MarkerSize',8,'MarkerFaceColor','b')
hold on
plot(time(6),zDot(vect(6)),'ko','MarkerSize',8,'MarkerFaceColor','k')
hold on
plot([time(1), time(6)], [p(1) p(1)], 'r-')
hold on
plot([time(6), time(11)], [p2(1) p2(1)], 'b-')
hold on
plot([time(6), time(6)], ylim, 'k--')
legend('Pre-Impact velocities', 'Post-impact velocities', 'Impact velocity','Pre-linear fit diff.','Post-linear fit diff.','Location','Southeast')
title({['Z velocity analysis for trial ', num2str(n)], '+- 5 data points'})
xlabel('Time (s)')
ylabel('Velocity (m/s)')
set(gca, 'FontSize', 12)


%All positions fitted!
figure 
tiledlayout(2,3)

%X position
nexttile
p3 = polyfit(time(1:6),x(vect(1:6))',1)
x1 = time(1):0.0001:time(6);
y3 = p3(1)*x1 + p3(2);
p4 = polyfit(time(6:11),x(vect(6:11))',1)
x2 = time(6):0.0001:time(11);
y4 = p4(1)*x2 + p4(2);
plot(time(1:5),x(vect(1:5)),'ro','MarkerSize',8,'MarkerFaceColor','r')
hold on
plot(time(7:11),x(vect(7:11)),'bo','MarkerSize',8,'MarkerFaceColor','b')
hold on
plot(time(6),x(vect(6)),'ko','MarkerSize',8,'MarkerFaceColor','k')
hold on
plot(x1,y3,'r')
hold on
plot(x2,y4,'b')
hold on
plot([time(6), time(6)], ylim, 'k--')
legend('Pre-Impact points', 'Post-impact points', 'Impact point','Pre-linear fit','Post-linear fit')
title({['X position analysis for trial ', num2str(n)], '+- 5 data points'})
xlabel('Time (s)')
ylabel('Position (m)')
set(gca, 'FontSize', 12)
xlim([time(1),time(11)])

%Z position
nexttile
p = polyfit(time(1:6),z(vect(1:6))',1)
x1 = time(1):0.0001:time(6);
y = p(1)*x1 + p(2);
p2 = polyfit(time(6:11),z(vect(6:11))',1)
x2 = time(6):0.0001:time(11);
y2 = p2(1)*x2 + p2(2);
plot(time(1:5),z(vect(1:5)),'ro','MarkerSize',8,'MarkerFaceColor','r')
hold on
plot(time(7:11),z(vect(7:11)),'bo','MarkerSize',8,'MarkerFaceColor','b')
hold on
plot(time(6),z(vect(6)),'ko','MarkerSize',8,'MarkerFaceColor','k')
hold on
plot(x1,y,'r')
hold on
plot(x2,y2,'b')
hold on
plot([time(6), time(6)], ylim, 'k--')
legend('Pre-Impact points', 'Post-impact points', 'Impact point','Pre-linear fit','Post-linear fit')
title({['Z position analysis for trial ', num2str(n)], '+- 5 data points'})
xlabel('Time (s)')
ylabel('Position (m)')
set(gca, 'FontSize', 12)
xlim([time(1),time(11)])

%Angular position
nexttile
p5 = polyfit(time(1:6),th(vect(1:6))',1)
x1 = time(1):0.0001:time(6);
y5 = p5(1)*x1 + p5(2);
p6 = polyfit(time(6:11),th(vect(6:11))',1)
x2 = time(6):0.0001:time(11);
y6 = p6(1)*x2 + p6(2);
plot(time(1:5),th(vect(1:5)),'ro','MarkerSize',8,'MarkerFaceColor','r')
hold on
plot(time(7:11),th(vect(7:11)),'bo','MarkerSize',8,'MarkerFaceColor','b')
hold on
plot(time(6),th(vect(6)),'ko','MarkerSize',8,'MarkerFaceColor','k')
hold on
plot(x1,y5,'r')
hold on
plot(x2,y6,'b')
hold on
plot([time(6), time(6)], ylim, 'k--')
legend('Pre-Impact points', 'Post-impact points', 'Impact point','Pre-linear fit','Post-linear fit')
title({['Angle analysis for trial ', num2str(n)], '+- 5 data points'})
xlabel('Time (s)')
ylabel('Angle (rads)')
set(gca, 'FontSize', 12)
xlim([time(1),time(11)])

%X velocity
nexttile
plot(time(1:5),xDot(vect(1:5)),'ro','MarkerSize',8,'MarkerFaceColor','r')
hold on
plot(time(7:11),xDot(vect(7:11)),'bo','MarkerSize',8,'MarkerFaceColor','b')
hold on
plot(time(6),xDot(vect(6)),'ko','MarkerSize',8,'MarkerFaceColor','k')
hold on
plot([time(1), time(6)], [p3(1) p3(1)], 'r-')
hold on
plot([time(6), time(11)], [p4(1) p4(1)], 'b-')
hold on
plot([time(6), time(6)], ylim, 'k--')
legend('Pre-Impact velocities', 'Post-impact velocities', 'Impact velocity','Pre-linear fit diff.','Post-linear fit diff.','Location','Southeast')
title({['X velocity analysis for trial ', num2str(n)], '+- 5 data points'})
xlabel('Time (s)')
ylabel('Velocity (m/s)')
set(gca, 'FontSize', 12)
xlim([time(1),time(11)])

%Z velocity
nexttile
plot(time(1:5),zDot(vect(1:5)),'ro','MarkerSize',8,'MarkerFaceColor','r')
hold on
plot(time(7:11),zDot(vect(7:11)),'bo','MarkerSize',8,'MarkerFaceColor','b')
hold on
plot(time(6),zDot(vect(6)),'ko','MarkerSize',8,'MarkerFaceColor','k')
hold on
plot([time(1), time(6)], [p(1) p(1)], 'r-')
hold on
plot([time(6), time(11)], [p2(1) p2(1)], 'b-')
hold on
plot([time(6), time(6)], ylim, 'k--')
legend('Pre-Impact velocities', 'Post-impact velocities', 'Impact velocity','Pre-linear fit diff.','Post-linear fit diff.','Location','Southeast')
title({['Z velocity analysis for trial ', num2str(n)], '+- 5 data points'})
xlabel('Time (s)')
ylabel('Velocity (m/s)')
set(gca, 'FontSize', 12)
xlim([time(1),time(11)])

%Angular Velocity
nexttile
plot(time(1:5),thDot(vect(1:5)),'ro','MarkerSize',8,'MarkerFaceColor','r')
hold on
plot(time(7:11),thDot(vect(7:11)),'bo','MarkerSize',8,'MarkerFaceColor','b')
hold on
plot(time(6),thDot(vect(6)),'ko','MarkerSize',8,'MarkerFaceColor','k')
hold on
plot([time(1), time(6)], [p5(1) p5(1)], 'r-')
hold on
plot([time(6), time(11)], [p6(1) p6(1)], 'b-')
hold on
plot([time(6), time(6)], ylim, 'k--')
legend('Pre-Impact velocities', 'Post-impact velocities', 'Impact velocity','Pre-linear fit diff.','Post-linear fit diff.','Location','Southeast')
title({['Angular velocity analysis for trial ', num2str(n)], '+- 5 data points'})
xlabel('Time (s)')
ylabel('Velocity (rads/s)')
set(gca, 'FontSize', 12)
xlim([time(1),time(11)])

end




