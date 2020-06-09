%% Compare Models 

%Try to debug our method of implementing the Wang Mason model with 
%Nima Fazeli's differantial equation approach
data = zeros(10,6);

%% Set up Variables
for i = 3:3
    n = i; %current trial for data being examined

    u = 0.071;  %mu, coefficiant of friction
    e = 0.568  %epsilon, coeeficiant of restitution 

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
    m = 36.4*1e-3; %Mass of ellipse [kg]

    %create Jacobian
    J = [d; n];
    %Nima's constants
    Minv = J*(m\J');
    M    = inv(Minv);

    a_e = 70/1000/2;
    b_e = 50/1000/2;
    %moment of inertia
    I_inertia = m*(b_e^2+a_e^2)/4;
    %general mass matrix
    Mass = [m, 0, 0;
            0, m, 0;
            0, 0, I_inertia];
    
    h = 1/250;
    ha = h*[0;-g;0];  


    %% Run Our Model
    out_juniors = APPoisson_juniors(Mass, n, d, pre(4:6), u ,e);

    %% Run Nima's Model
    [out_nima, z] = APPoisson_nima(Mass, n', d', pre(4:6)', ha, u, e);

    %% Compare Results
    disp("Actual Velocity")
    disp(post(4:5))
    disp("Our Model")
    disp(out_juniors)
    disp("Nima's Model")
    disp(out_nima)

    data(i, 1) = post(4);
    data(i, 2) = out_juniors(1);
    data(i, 3) = out_nima(1)';
    data(i, 4) = post(5);
    data(i, 5) = out_juniors(2);
    data(i, 6) = out_nima(2)';
    data(i, 7) = post(6);
    data(i, 8) = out_juniors(3);
    data(i, 9) = out_nima(3);
    
end