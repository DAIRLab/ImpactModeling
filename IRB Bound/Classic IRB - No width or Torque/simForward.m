trial = 585;

pos = bounce_array(trial).states(1:3)';
vel = bounce_array(trial).states(4:6)';

next = true; 
dt = 1/250/100;
totalTime = 0;

%% Simulate ellipse forward in time until part of the ellipse contacts ground
J_0 = getJacobian(pos);
while next == true
    %Update COM State
    x_com = pos(1) + vel(1) * dt; %x = x_0 + v_x * t
    y_com = pos(2) + vel(2) * dt + 4.9 * dt^2; % y = y_0 + v_y * t + 1/2 g * t^2
    th_com = pos(3) + vel(3) * dt; %th = th_0 + w * t
    
    v_ycom = vel(2) + 9.8 * dt; %v_y = v_y0 + g * t
    
    %update position\velocity vectors
    pos = [x_com; y_com; th_com];
    vel = [vel(1); v_ycom; vel(3)];
    
    %find points on ellipse
    [xvec,yvec] = ellipse_visual(pos);
    
    %check if contact occurs
    if min(yvec) <= 0
        next = false;
    end
    
    totalTime = totalTime + dt;
    
end 
J_f = getJacobian(pos);
disp(pos)
disp(vel)
disp("Original J")
disp(J_0)
disp("Final J")
disp(J_f)