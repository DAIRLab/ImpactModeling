%% Finding Contact Jacobian

% Function: Finds the normal and tangential components of the Contact
%           Jacobian by including both position and angular components in
%           calculations.
%           **Only finds first impact for now
% Input: tr - trial number
%        impact - impact to analyze
% Output: d - tangential 
%         n - normal

function [n, d] = findNDupdated(tr, impact)

%load data
traj = readtable("traj_" + tr + ".csv");    
%convert from table to matrix
traj = traj{:, :};


side = 0.06/2; %meters
dt = 1/10000000; %seconds (must be less than 1/250 because the framerate is 250Hz)
t = 0; %start time at 0
g = -9.8; %m/s^2
x = traj(:,1); %Column 1: x position
y = traj(:,2);  %Column 2: y position
th = traj(:,3);  %Column 3: theta
vx = traj(:,4);  %Column 4: x velocity --> constant
vy = traj(:,5);  %Column 5: y velocity --> increasing due to gravity
w = traj(:,6);  %Column 6: angular velocity --> assume constant ?
% tol = 0.00088; 
tol = 0.00005;
buffer = 2; %how many frames to go back


imp = zeros(1, 6); %[x, y, th, vx, vy, w]

% initialize n and d 
n = [0, 1, 0]; %eventually [0, 1, xcom - xcontact]
d = [1, 0, 0]; %eventually [1, 0, ycom - ycontact]

% 1. extrapolate the impact 
% 2. use the preimpact data
% 3. assume x velocity to be constant and update x (x = vx*t + x)
% 4. increase y velocity according to gravity (y = -0.5*g*t^2 + vy*t +
%    y)
% 5. assume angular velocity to be constant and update th (th = w*t + th)
% 6. update the square position by time steps until the lowest point is
%    zero/within a tolerance (update time --> t = t+dt)
% 7. calculate d and n by subtracting the com position from contact
%    position

% EXTRAPOLATE IMPACT
  count = 1;
  impactCount = 0;
  while(sum(imp)==0) %while impact vector has not been modified
      s0 = sign(vy(count));
      if (s0 == -1) && (sign(vy(count+1)) == 1) %checking for impact
          impactCount = impactCount + 1;
          if impact == impactCount  %checking for correct impact
              % use pre-impact data  --> FOR NOW: try going back a frame
              imp(1) = x(count - buffer); 
              imp(2) = y(count - buffer);
              imp(3) = th(count - buffer);
              imp(4) = vx(count - buffer);
              imp(5) = vy(count - buffer);
              imp(6) = w(count - buffer);
          end
      else
         %move on to next row
         count = count + 1;
      end
  end

  xPos = imp(1); %com
  yPos = imp(2); %com
  theta = imp(3);
  xVel = imp(4);
  yVel = imp(5);
  omega = imp(6);
  
  %START andy's method of calculating corners

     c1 = xPos;
     c2 = yPos;
     ax = side;
     ay = side;
     bx = side;
     by = -side;
     cx = -side;
     cy = -side;
     dx = -side;
     dy = side;

     a = [ax;ay];
     b = [bx;by];
     c = [cx;cy];
     f = [dx;dy];
     
     cornersX = [a(1),b(1),c(1),f(1),a(1)];
     cornersY = [a(2),b(2),c(2),f(2),a(2)];

     R = [cos(theta),sin(theta);-sin(theta),cos(theta)];

     a_2 = R*a;
     b_2 = R*b;
     c_2 = R*c;
     f_2 = R*f;

     cornersX = [a_2(1),b_2(1),c_2(1),f_2(1),a_2(1)]+c1;
     cornersY = [a_2(2),b_2(2),c_2(2),f_2(2),a_2(2)]+c2;
     
     ycontact = min(cornersY);
     xcontact = cornersX(find(ycontact==cornersY));
     
  % CALCULATE EXACT IMPACT STATE
  while((abs(ycontact) > tol) && (t <= buffer/250)) %terminate if exceeds buffer frames
     xPos = xVel*t + xPos;
     yPos = 0.5*g*(t^2) + yVel*t + yPos;
     yVel = g*t + yVel;
     theta = omega*t + theta; 
      
     c1 = xPos;
     c2 = yPos;
     ax = side;
     ay = side;
     bx = side;
     by = -side;
     cx = -side;
     cy = -side;
     dx = -side;
     dy = side;

     a = [ax;ay];
     b = [bx;by];
     c = [cx;cy];
     f = [dx;dy];

     cornersX = [a(1),b(1),c(1),f(1),a(1)];
     cornersY = [a(2),b(2),c(2),f(2),a(2)];

     R = [cos(theta),sin(theta);-sin(theta),cos(theta)];

     a_2 = R*a;
     b_2 = R*b;
     c_2 = R*c;
     f_2 = R*f;

     cornersX = [a_2(1),b_2(1),c_2(1),f_2(1),a_2(1)]+c1;
     cornersY = [a_2(2),b_2(2),c_2(2),f_2(2),a_2(2)]+c2;

     ycontact = min(cornersY);
     xcontact = cornersX(find(ycontact==cornersY));
     xcontact = xcontact(1);
     %update t
     t = dt + t;
  end

  if t > buffer/250 %if solution was not found, output zeros
    n = zeros(1,3);
    d = n;
  else
    d(3) = yPos - ycontact;
    n(3) = xcontact - xPos;
  end
end
