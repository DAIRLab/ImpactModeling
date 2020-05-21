
%{
 This code reads in the .csv data file and creates plots of the trajectory
 of the square (both position vs time and x vs y). Additionally, the 
 points of impact are marked by red squares on the plots. The radius of
 gyration is calculated based on the moment of inertia of a square and the 
 given dimension of the square.
%}

s = 0.06; %m (side length of square)

rad_gyr = sqrt(s^2 / 6); %radius of gyration of square - m

disp("radius of gyration is " + Izz + " meters");

traj1 = csvread('traj_1.csv', 1, 0);

x = traj1(:,1); %m
y = traj1(:,2); %m
theta = traj1(:,3); %radians ?
vx = traj1(:,4); %m/s ?
vy = traj1(:,5); %m/s ?
vtheta = traj1(:,6); %radians/s ?

tstep = 1/250; %seconds
t = zeros(200,1);

impactX = zeros();
impactY = zeros();
impactT = zeros();

for i = 1:200
    %create time vector
    t(i,:) = (i-1)*tstep;
    
    %set intermediate current value of y velocity to compare sign of next
    curr = vy(i);
    if (i < 200)
        next = vy(i+1);
    end
    
    %determine when there is an impact/"bounce"
    if (sign(curr) == -1 * sign(next))
        impactX(end+1) = x(i,:);
        impactY(end+1) = y(i,:);
        impactT(end+1) = t(i,:);
    end
    
end

%updating to delete the origin data point
impactX = impactX(2:end);
impactY = impactY(2:end);
impactT = impactT(2:end);


figure(1)
plot(t, y)
hold on 
plot(impactT, impactY, 's', 'MarkerSize', 10)
title('Y Position over Time')
xlabel('Time (s)')
ylabel('Y Position (m)')

figure(2)
plot(x,y)
hold on
plot(impactX, impactY, 's', 'MarkerSize' , 10)
title('Trajectory')
xlabel('X Position (m)')
ylabel('Y Position (m)')


