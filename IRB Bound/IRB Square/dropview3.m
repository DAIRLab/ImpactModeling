%convert from table to matrix
height = 0.06/2;
width = 0.06/2;

% Visualize the Drop angle
traj = readmatrix(['/Users/andyeske/Desktop/Simulations/dice-data-processed/traj_' num2str(6) '.csv']);   

x = traj(:,1);
%Column 2: z position
z = traj(:,2);
%Column 3: theta
th = traj(:,3);
%Column 4: xDot
xDot = traj(:,4);
%Column 5: zDot
zDot = traj(:,5);
%Column 6: thetaDot
thDot = traj(:,6);

% Visualize the Drop angle2
traj2 = readmatrix(['/Users/andyeske/Desktop/Simulations/simulated-trajectories/trajSim_' num2str(6) '.csv']);   

x2 = traj2(:,1);
%Column 2: z position
z2 = traj2(:,2);
%Column 3: theta
th2 = traj2(:,3);
%Column 4: xDot
xDot2 = traj2(:,4);
%Column 5: zDot
zDot2 = traj2(:,5);
%Column 6: thetaDot
thDot2 = traj2(:,6);

figure
tiledlayout(1,2)

nexttile

for i = 1:(length(traj)-30)
     hold on
     plot(-x(i),z(i),'r*')
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
     ylim([-0.1 0.25])
     xlim([0.2 0.75])
     plot(vx,vy, 'k')
     axis equal
     axis off
     pause(1/250);

end

nexttile

for i = 1:(length(traj2)-30)
     hold on
     plot(-x2(i),z2(i),'r*')
     the = th2(i);
     c1_2 = -x2(i);
     c2_2 = z2(i);
     ax_2 = 0.03;
     ay_2 = 0.03;
     bx_2 = 0.03;
     by_2 = - 0.03;
     cx_2 = - 0.03;
     cy_2 = - 0.03;
     dx_2 = - 0.03;
     dy_2 = 0.03;

     a2 = [ax_2;ay_2];
     b2 = [bx_2;by_2];
     c2 = [cx_2;cy_2];
     d2 = [dx_2;dy_2];
     
     vx2 = [a2(1),b2(1),c2(1),d2(1),a2(1)];
     vy2 = [a2(2),b2(2),c2(2),d2(2),a2(2)];

     R2 = [cos(the),sin(the);-sin(the),cos(the)];

     a_2 = R2*a2;
     b_2 = R2*b2;
     c_2 = R2*c2;
     d_2 = R2*d2;

     vx2 = [a_2(1),b_2(1),c_2(1),d_2(1),a_2(1)]+c1_2;
     vy2 = [a_2(2),b_2(2),c_2(2),d_2(2),a_2(2)]+c2_2;
     ylim([-0.1 0.25])
     xlim([0.2 0.75])
     plot(vx2,vy2, 'k')
     axis equal
     axis off
     pause(1/250);

end
