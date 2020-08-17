%% Finding Contact Jacobian

% Function: Finds the normal and tangential components of the Contact
%           Jacobian by including both position and angular components in
%           calculations.
%           **Only finds first impact for now
% Input: tr - trial number
% Output: d - tangential 
%         n - normal

function [n, d] = findND(tr)

%load data
traj = readtable("traj_" + tr + ".csv");    
%convert from table to matrix
traj = traj{:, :};


height = 0.06/2; %meters
width = 0.06/2; %meters
dt = 1/7500; %seconds (must be less than 1/250 because the framerate is 250Hz)
t = 0; %start time at 0
g = -9.8; %m/s^2
x = traj(:,1); %Column 1: x position
y = traj(:,2);  %Column 2: y position
th = traj(:,3);  %Column 3: theta
vx = traj(:,4);  %Column 4: x velocity --> constant
vy = traj(:,5);  %Column 5: y velocity --> increasing due to gravity
w = traj(:,6);  %Column 6: angular velocity --> assume constant ?
tol = 0.001; 

impact1 = zeros(1, 6); %[x, y, th, vx, vy, w]

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
  while(sum(impact1)==0) %while impact vector has not been modified
      s0 = sign(vy(count));
      if (s0 == -1) && (sign(vy(count+1)) == 1)
          % use pre-impact data
          impact1(1) = x(count); 
          impact1(2) = y(count);
          impact1(3) = th(count);
          impact1(4) = vx(count);
          impact1(5) = vy(count);
          impact1(6) = w(count);
      else
         %move on to next row
         count = count + 1;
      end
  end

  xPos = impact1(1); %com
  yPos = impact1(2); %com
  theta = impact1(3);
  xVel = impact1(4);
  yVel = impact1(5);
  omega = impact1(6);
  
  R = [cos(theta), -sin(theta), sin(theta), cos(theta)];
  cX = [-width, width, width, -width];
  cY = [-height, -height, height, height];
  T(1,:)=R.*cX;
  T(2,:)=R.*(cY);

  
  %initialize corners
  xll = xPos + T(1,1); %lower left
  xlr = xPos + T(1,2); %lower right
  xur = xPos + T(1,3); %upper left
  xul = xPos + T(1,4); %upper left
  yll = yPos + T(2,1); 
  ylr = yPos + T(2,2);
  yur = yPos + T(2,3);
  yul = yPos + T(2,4);
  
  cornersX = [xll, xlr, xur, xul];
  cornersY = [yll, ylr, yur, yul];
  
  ycontact = min(cornersY);
  xcontact = cornersX(find(ycontact==cornersY)); %corresponding xcontact

  % CALCULATE EXACT IMPACT STATE
  while(abs(ycontact) > tol)
      xPos = xVel*t + xPos;
      yPos = 0.5*g*(t^2) + yVel*t + yPos;
      yVel = g*t + yVel;
      theta = omega*t + theta; 
      
      R = [cos(theta), -sin(theta), sin(theta), cos(theta)];
      T(1,:)=R.*cX;
      T(2,:)=R.*(cY);
      
      xll = xPos + T(1,1);
      xlr = xPos + T(1,2);
      xur = xPos + T(1,3);
      xul = xPos + T(1,4);
      yll = yPos + T(2,1);
      ylr = yPos + T(2,2);
      yur = yPos + T(2,3);
      yul = yPos + T(2,4);

      cornersX = [xll, xlr, xur, xul];
      cornersY = [yll, ylr, yur, yul];

      ycontact = min(cornersY);
      xcontact = cornersX(find(ycontact==cornersY)); %corresponding xcontact
  
      %update t
      t = dt + t;
  end

  d(3) = yPos - ycontact;
  n(3) = xPos - xcontact;
  
end
