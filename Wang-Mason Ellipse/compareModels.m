%% Compare Models 

%Try to debug our method of implementing the Wang Mason model with 
%Nima Fazeli's differantial equation approach
data = zeros(10,6);

%% Set up Variables
for i = 1:30
    n = i; %current trial for data being examined

    u = 0.057;  %mu, coefficiant of friction
    e = 0.612;  %epsilon, coeeficiant of restitution 

    load('ellipse_uniform.mat'); %load in ellipse collision data
    %pre and post impact data
    [pre] = bounce_array(n).states(1:6); 
    [post] = bounce_array(n).states(7:12);
    d = bounce_array(n).d;
    n = bounce_array(n).n;

    %assign pre impact positions for object 1
    th = pre(1,3);
    y1 = pre(1,2); 


    g = 9.81; %gravitational constant
    mass = 36.4*1e-3; %Mass of ellipse [kg]


    %create Jacobian
    J = [d; n];
    %Nima's constants
    Minv = J*(mass\J');
    M    = inv(Minv);

    a_e = 70/1000/2;
    b_e = 50/1000/2;
    I_inertia = mass*(b_e^2+a_e^2)/5;
    Mass = [mass,0,0;
            0,mass,0;
            0,0,I_inertia];
    h = 1/250;
    ha = h*[0;-g;0];
    
    %Find rho to use in equations for B1 -> B3
    rho = sqrt(I_inertia/mass); % sqrt(I/m)    

    %%  Solve for constants
    %find the sign of x1, ie whether it is on the left or right of the contact point.
    if th > 0 %turning anticlockwise
        r = rem(th,pi)/pi;
        if r > 0 && r < 0.25 || r > 0.5 && r < 0.75 %angle between 0 to 90 degrees and 180-270
            signx1 = 1; %x1 is to the right of contact point
        elseif r > 0.25 && r < 0.5 || r > 0.75 %angle between 90 and 180 and 270-360
            signx1 = -1; %x1 will be on the left
        else %when it is right on top
            signx1 = 1;
        end
    elseif th < 0 %turning clockwise
        r = rem(abs(th),pi)/pi;
        if r > 0 && r < 0.25 || r > 0.5 && r < 0.75 %angle between 0 to 90 degrees and 180-270
            signx1 = -1; %x1 is to the left of contact point
        elseif r > 0.25 && r < 0.5 || r > 0.75 %angle between 90 and 180 and 270-360
            signx1 = 1; %x1 will be on the right
        else %when it is right on top
            signx1 = 1;
        end
    elseif th == 0 %doesn't even rotate
        signx1 = 1;
    else
        signx1 = 1; 
    end

    %use a tilted ellipse equation to find x1 (x distance between the COM and the contact point)
    a = 0.035;
    b = 0.025;
    k = 0;
    h = 0;
    y = -y1;

%     x1 = double((a*b*(b^2*cos(th)^2 - k^2*cos(th)^4 + a^2*sin(th)^2 ...
%      - y^2*cos(th)^4 - k^2*sin(th)^4 - y^2*sin(th)^4 + 2*k*y*cos(th)^4 + ...
%      2*k*y*sin(th)^4 - 2*k^2*cos(th)^2*sin(th)^2 - ...
%      2*y^2*cos(th)^2*sin(th)^2  + 4*k*y*cos(th)^2*sin(th)^2)^(1/2) +...
%      b^2*h*cos(th)^2 + a^2*h*sin(th)^2 - a^2*k*cos(th)*sin(th) + ...
%      b^2*k*cos(th)*sin(th) + a^2*y*cos(th)*sin(th) -...
%      b^2*y*cos(th)*sin(th))/(a^2*sin(th)^2 + b^2*cos(th)^2));
% 
%     x1 = signx1*abs(x1(1)); 

    %USING POSA'S SUGESTION
    x1 = n(3);
    y1 = d(3);
    
    B1 = 1/mass + y1^2/(mass*rho^2); %(19)
    B2 = 1/mass + x1^2/(mass*rho^2); %(20)
    B3 = x1 * -y1/mass/rho^2;    %(21)

    V_c = J*[pre(1,4);pre(1,5);pre(1,6)];
    S_0 = V_c(1); %(22)
    C_0 = V_c(2); %(23)


    %% Run Our Model
    out_juniors = wang_juniors(mass, S_0, C_0, [B1; B2; B3], pre(4:6), u ,e);


    %% Run Nima's Model
    [out_nima, z] = wang_nima(Mass, n', d', pre(4:6)', ha, u, e);


    %% Compare Results
    disp("Actual Velocity")
    disp(post(4:5))
    disp("Our Model")
    disp(out_juniors)
    disp("Nima's Model")
    disp(out_nima(1:2))
    % norm(out_juniors - post(4:5))
    % norm(out_nima - post(4:5))
    data(i, 1) = post(4);
    data(i, 2) = out_juniors(1);
    data(i, 3) = out_nima(1)';
    data(i, 4) = post(5);
    data(i, 5) = out_juniors(2);
    data(i, 6) = out_nima(2)';
    
end