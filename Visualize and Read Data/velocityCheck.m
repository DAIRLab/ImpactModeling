%Manually checking accuracy of velocity for square drop data set
%1. Graph positions around impact and calculate velocities from slopes of 
%fit curves
%2. Compare to velocity in actual data set

%% Impact 6
traj = readtable("traj_6.csv");    
%convert from table to matrix
traj = traj{:, :};
rate = 1/250; %capture rate
%impact occurs between line 15-16
impact = 14;
%plot 10 frames before and after
range = impact-10:impact+10;
time = (impact-10)*rate:rate:(impact+10)*rate;  %time vector
position_x = traj(range, 1);                    %position vector
velocities_x = traj(range, 4);                  %velcoity vector
figure()
subplot(2,1,1)
hold on
%look at x velocities
plot(time, position_x, 'k.')
%plot impact points
plot(time(11:12), position_x(11:12, 1), 'k*')


%pre-impact trend
p1 = polyfit(time(1:11)', position_x(1:11), 1); 
plot(time(1:11), p1(1)*time(1:11)+p1(2),'b-');
v_pre = p1(1);
leg1 = ['Pre Impact Trend, m = '  num2str(v_pre)];

%post-impact trend
p2 = polyfit(time(12:end)', position_x(12:end), 1); 
plot(time(12:end), p2(1)*time(12:end)+p2(2), 'r-')
v_post = p2(1);
leg2 = ['Post Impact Trend, m = '  num2str(v_post)];
xlabel("Time [s]")
ylabel("X Position [m]")
title("X Position")
legend('X Positions Observed', 'Pre and Post Impact X Positions', leg1, leg2)

%%
subplot(2,1,2)
hold on
plot(time, velocities_x, 'k.');
plot(time(11:12), velocities_x(11:12, 1), 'k*')
plot(time(1:11), ones(11,1)*v_pre, 'b-');
plot(time(12:end), ones(10,1)*v_post, 'r-');
xlabel("Time [s]")
ylabel("X Velocity [m/s]")
title("X Velocity")
legend('X Velocities Observed', 'X Impact Velocities', 'Predicted Pre-Impact X Velocity', ...
    'Predicted Post-Impact X Velocity', 'Location', 'Southeast');

%%
position_y = traj(range, 2);                    %position vector
velocities_y = traj(range, 5);                  %velcoity vector
figure()
subplot(2,1,1)
hold on
%look at x velocities
plot(time, position_y, 'k.')
%plot impact points
plot(time(11:12), position_y(11:12, 1), 'k*')


%pre-impact trend
p1 = polyfit(time(1:11)', position_y(1:11), 1); 
plot(time(1:11), p1(1)*time(1:11)+p1(2),'b-');
v_pre = p1(1);
leg1 = ['Pre Impact Trend, m = '  num2str(v_pre)];

%post-impact trend
p2 = polyfit(time(12:end)', position_y(12:end), 1); 
plot(time(12:end), p2(1)*time(12:end)+p2(2), 'r-')
v_post = p2(1);
leg2 = ['Post Impact Trend, m = '  num2str(v_post)];
xlabel("Time [s]")
ylabel("Y Position [m]")
title("Y Position")
legend('Y Positions Observed', 'Pre and Post Impact Y Positions', leg1, leg2)

%%
subplot(2,1,2)
hold on
plot(time, velocities_y, 'k.');
plot(time(11:12), velocities_y(11:12, 1), 'k*')
plot(time(1:11), ones(11,1)*v_pre, 'b-');
plot(time(12:end), ones(10,1)*v_post, 'r-');
xlabel("Time [s]")
ylabel("Y Velocity [m/s]")
title("Y Velocity")
legend('Y Velocities Observed', 'Y Impact Velocities', 'Predicted Pre-Impact Y Velocity', ...
    'Predicted Post-Impact Y Velocity', 'Location', 'Southeast');
%%
position_theta = traj(range, 3);                    %position vector
velocities_theta = traj(range, 6);                  %velcoity vector
figure()
subplot(2,1,1)
hold on
%look at x velocities
plot(time, position_theta, 'k.')
%plot impact points
plot(time(11:12), position_theta(11:12, 1), 'k*')


%pre-impact trend
p1 = polyfit(time(1:11)', position_theta(1:11), 1); 
plot(time(1:11), p1(1)*time(1:11)+p1(2),'b-');
v_pre = p1(1);
leg1 = ['Pre Impact Trend, m = '  num2str(v_pre)];

%post-impact trend
p2 = polyfit(time(12:end)', position_theta(12:end), 1); 
plot(time(12:end), p2(1)*time(12:end)+p2(2), 'r-')
v_post = p2(1);
leg2 = ['Post Impact Trend, m = '  num2str(v_post)];
xlabel("Time [s]")
ylabel("Theta [rad]")
title("Theta")
legend('Theta Observed', 'Pre and Post Impact Theta', leg1, leg2)

%%
subplot(2,1,2)
hold on
plot(time, velocities_theta, 'k.');
plot(time(11:12), velocities_theta(11:12, 1), 'k*')
plot(time(1:11), ones(11,1)*v_pre, 'b-');
plot(time(12:end), ones(10,1)*v_post, 'r-');
xlabel("Time [s]")
ylabel("Theta Dot [rad/s]")
title("Theta Dot")
legend('Theta Observed', 'Theta Velocities', 'Predicted Pre-Impact Theta Dot', ...
    'Predicted Post-Impact Theta Dot', 'Location', 'Southeast');

