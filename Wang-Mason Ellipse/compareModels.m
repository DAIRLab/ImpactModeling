%% Compare Models 

%Try to debug our method of implementing the Wang Mason model with 
%Nima Fazeli's differantial equation approach

%% Set up Variables
n = 5; %current trial for data being examined

u = 0.057; %mu, coefficiant of friction
e = 0.612;  %epsilon, coeeficiant of restitution 

load('ellipse_uniform.mat'); %load in ellipse collision data
%pre and post impact data
[pre] = bounce_array(n).states(1:6); 
[post] = bounce_array(n).states(7:12);
d = bounce_array(n).d;
n = bounce_array(n).n;

%assign pre impact velocities for object 1
th = pre(1,3);
y1 = pre(1,2); 

%Find rho
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
rho = 0.5 * sqrt(a0^2 + b0^2); % sqrt(I/m)

g = 9.81; %gravitational constant
mass = 1; %Mass of ellipse [kg]
h = 1/250;
ha = h*[0;-g;0];

%create Jacobian
J = [d; n];

% Minv = J*(Mass\J');
% M    = inv(Minv);

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

x1 = double((a*b*(b^2*cos(th)^2 - k^2*cos(th)^4 + a^2*sin(th)^2 ...
 - y^2*cos(th)^4 - k^2*sin(th)^4 - y^2*sin(th)^4 + 2*k*y*cos(th)^4 + ...
 2*k*y*sin(th)^4 - 2*k^2*cos(th)^2*sin(th)^2 - ...
 2*y^2*cos(th)^2*sin(th)^2  + 4*k*y*cos(th)^2*sin(th)^2)^(1/2) +...
 b^2*h*cos(th)^2 + a^2*h*sin(th)^2 - a^2*k*cos(th)*sin(th) + ...
 b^2*k*cos(th)*sin(th) + a^2*y*cos(th)*sin(th) -...
 b^2*y*cos(th)*sin(th))/(a^2*sin(th)^2 + b^2*cos(th)^2));

x1 = signx1*abs(x1(1)); 

B1 = 1/mass + y1^2/(mass*rho^2); %(19)
B2 = 1/mass + x1^2/(mass*rho^2); %(20)
B3 = x1 * y1/mass/rho^2;    %(21)

V_c = J*[pre(1,4);pre(1,5);pre(1,6)];
S_0 = V_c(1); %(22)
C_0 = V_c(2); %(23)


%% Run Our Model
out_juniors = wang_juniors(mass, S_0, C_0, [B1; B2; B3], pre(4:6), u ,e);


%% Run Nima's Model
[out_nima, z] = wang_nima(mass, n', d', pre(4:6)', ha, u, e);


%% Compare Results
