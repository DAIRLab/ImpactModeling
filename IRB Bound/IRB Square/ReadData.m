function [Pre, Post] = ReadData(n)

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
%plot(-x,z)

for jj = 1:l
     %hold on
     %pause(0.05)
     %plot(-x(jj),z(jj),'o')
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
%     
     %plot(vx,vy)

     R = [cos(the),sin(the);-sin(the),cos(the)];

     a2 = R*a;
     b2 = R*b;
     c2 = R*c;
     d2 = R*d;

     vx = [a2(1),b2(1),c2(1),d2(1),a2(1)]+cc1;
     vy = [a2(2),b2(2),c2(2),d2(2),a2(2)]+cc2;
    
     %plot(vx,vy)
     %ylim([-0.1 0.25])
     %xlim([0.2 0.75])
     %axis equal
end

% for jj = 1:10
%      hold on
%      pause(0.1)
% 
% the = th(jj);
% vx = [0, cos(the)];
% vy = [0, sin(the)];
% 
% plot(vx,vy)
% 
% end

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

Pre = [impR', x(impR), z(impR), th(impR), xDot(impR), zDot(impR), thDot(impR)];
Post = [impR2', x(impR2), z(impR2), th(impR2), xDot(impR2), zDot(impR2), thDot(impR2)];

%plotting for all of the elements without assuming any time span
% t = 0:1:l-1;
% plot(t,z)
% title("Tracking changes in z position");
% ylabel("z position")
% xlabel("time");
% hold on
% plot(impR,z(impR),'o')

% the = pi/6;
% c1 = 0.3
% c2 = 0.3
% ax = c1 + 0.03;
% ay = c2 + 0.03;
% bx = c1 + 0.03;
% by = c2 - 0.03;
% cx = c1 - 0.03;
% cy = c2 - 0.03;
% dx = c1 - 0.03;
% dy = c2 + 0.03;
% 
% a = [ax;ay];
% b = [bx;by];
% c = [cx;cy];
% d = [dx;dy];
% 
% xlim([0 .4])
% ylim([0 .4])
% 
% R = [cos(the),sin(the);-sin(the),cos(the)];
% 
% a2 = R*a;
% b2 = R*b;
% c2 = R*c;
% d2 = R*d;
% 
% vx = [a2(1),b2(1),c2(1),d2(1),a2(1)];
% vy = [a2(2),b2(2),c2(2),d2(2),a2(2)];
% plot(vx,vy)
end











