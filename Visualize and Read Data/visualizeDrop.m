    %% Visualize the Drop angle
    %bigErr = [260,2; 268,2];
    trials = 191;
    figure();

for j = 14%badTrials(end)
    pause(1)
    clf;
    traj = readtable("traj_"+num2str(j)+".csv");    
    %convert from table to matrix
    traj = traj{:, :};
    height = 0.06/2;
    width = 0.06/2;

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
    for i = 1:(length(traj)-80)
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
%          text = ["Traj:" num2str(dataSet(j).trial) " Impact:" num2str(dataSet(j).impact)];
%          title(text);
         pause(10/250);

    end
    pause(5);
end