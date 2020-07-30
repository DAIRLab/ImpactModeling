%% Visualize simulation and drop together

figure()

for j = 1:20
    pause(1)
    %load trajectory csv files
    traj_real = readtable("traj_"+num2str(j) +".csv");   
    traj_sim = readtable("trajSim_"+num2str(j) +".csv"); 
    %convert from table to matrix
    traj_real = traj_real{:, :};
    traj_sim = traj_sim{:, :};
    
    %info about square
    height = 0.06/2;
    width = 0.06/2;
    
    %Column 2: x position
    x_real = traj_real(:,1);
    x_sim = traj_sim(:,1);
    %Column 2: z position
    z_real = traj_real(:,2);
    z_sim = traj_sim(:,2);
    %Column 3: theta
    th_real = traj_real(:,3);
    th_sim = traj_sim(:,3);
    %Column 4: xDot
    xDot_real = traj_real(:,4);
    xDot_sim = traj_sim(:,4);
    %Column 5: zDot
    zDot_real = traj_real(:,5);
    zDot_sim = traj_sim(:,5);
    %Column 6: thetaDot
    thDot_real = traj_real(:,6);
    thDot_sim = traj_sim(:,6);
    clf; 
    sgtitle(["Trial Number " num2str(j)])
    for i = 1:(length(traj_real)-80)
        %----------REAL TRAJECTORY-----
         subplot(2,1,1);
         hold on
         plot(-x_real(i),z_real(i),'r*')
         the = th_real(i);
         c1 = -x_real(i);
         c2 = z_real(i);
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
         text = ["Observed Trajectory"];
         title(text);
         
         %------SIMULATED TRAJECTORY-------
         subplot(2,1,2);
         hold on
         plot(-x_sim(i),z_sim(i),'g*')
         the = th_sim(i);
         c1 = -x_sim(i);
         c2 = z_sim(i);
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
         text = ["Simulated Trajectory"];
         title(text);
         
         pause(1/250);
            
    end
    pause(1);
end