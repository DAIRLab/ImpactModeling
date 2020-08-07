clear; 
load('ellipse_uniform.mat'); %load in ellipse collision data

numTrials = 200; %Number of Trials
muRange = linspace(0,1);
epsilonRange = linspace(0,1);

g = 9.81; %gravitational constant
m = 36.4*1e-3; %Mass of ellipse [kg]
%moment of inertia
a_e = 70/1000/2;
b_e = 50/1000/2;
I_inertia = m*(b_e^2+a_e^2)/4;
%general mass matrix
Mass = [m, 0, 0;
        0, m, 0;
        0, 0, I_inertia];

h = 1/250;
ha = h*[0;-g;0];  

error = zeros(100);
totalError = zeros(100);
count = 0;

for t = 1:numTrials
    trial = t;
    clear bestA;
    clear bestB;
    if sum(bounce_array(trial).flags)<1
    count = count + 1;
    %Get pre and post impact data from current trial
    pre = bounce_array(trial).states(4:6);
    post = bounce_array(trial).states(10:12);
    
    %get normal and tangetnial Jacobians from dataset
    d = (bounce_array(trial).d);   %tangential
    n = (bounce_array(trial).n);   %normal
    
    %create Jacobian
    J = [d; n];
    
        
    for a = 1:100
        mu = muRange(a);
        for b = 1:100
            epsilon = epsilonRange(b);
            out_juniors = APPoisson_juniors(Mass, n, d, pre, mu, epsilon);
            %[out_nima, z] = APPoisson_nima(Mass, n', d', pre', ha, mu, epsilon);
            error(a,b) = norm(post - out_juniors)/norm(post);
            %error(a,b) = norm(post - out_nima)/norm(post);
        end
    end
    minimum = min(min(error));
    [bestA, bestB]=find(error==minimum);
    bestMu(count) = muRange(min(bestA));
    bestEpsilon(count) = epsilonRange(min(bestB));
    totalError = totalError + error;
    end
end

disp(mean(bestMu))
disp(mean(bestEpsilon))

totalError = totalError / count; 
%%
contourf(epsilonRange, muRange, totalError', 30);
colorbar;
xlabel('Epsilon')
ylabel('Mu')
%title('Error for Mu and Epsilon over 20 Random Trials')