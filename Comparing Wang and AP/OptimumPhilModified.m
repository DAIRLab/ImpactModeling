% Find Best Mu and Epsilon using AP Poisson

% load in ellipse data 
load('ellipse_uniform.mat');

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
M = [m1,0,0;0,m1,0;0,0,I1];

% Set up interval
iter = 50; %how many iterations of mu and epsilon we would like
sample = linspace(0.01, 0.99, iter); %set up sample vector
trialError = zeros(iter);
totalError = zeros(iter);

% Run Simulation
%loop over multiple trials
numTrials = 500; %Number of Trials
bev = zeros(1,numTrials);
bmv = zeros(1,numTrials);
v = randi([1 2000],1,numTrials);

for tr = 1:numTrials
    trial = v(tr);
    %Get pre and post impact data from current trial
    pre = bounce_array(trial).states(4:6);
    post = bounce_array(trial).states(10:12);
    
    %get normal and tangetnial Jacobians from dataset
    d = (bounce_array(trial).d);   %tangential
    n = (bounce_array(trial).n);   %normal
    
    J = [d;n];
    Minv = J*(M\J');
    Ma = inv(Minv);
    
    %loop over all mu's
    for i = 1:iter
        %assign epsilon to value corresponding with loop
        epsilon = sample(i);

        %loop over all epsilons
        for j = 1:iter
            %assign mu to value corresponding with loop
            mu = sample(j);

            %Run AP Poisson Model given mu and epsilon
             v_calc = APPoisson_juniors(M, n, d, pre, mu, epsilon);
             
             Pt = (Ma*(J*(post-pre)'));
             p_hat = (Ma*(J*(v_calc-pre')));
             
             error = norm(Pt-p_hat);
             error = error/norm(Pt);
            
             trialError(i, j) = error;            
            
        end
    end
    %add trial error to total error
    totalError = totalError +  trialError;
    [be , bm] = find(trialError == min(min(trialError)));
    bev(tr) = min(be);
    bmv(tr) = min(bm);
end

% Post Process Results
%Take the average of totalError
averageError = totalError / numTrials;
minError = min(min(averageError));  %get the minimum error
[a , b] = find(averageError == minError); %find its index
%convert the index to mu and epsilon values
bestEp = sample(a);
bestMu = sample(b);

%Finding average mu and ep by taking the average of the best of each trial
aa = mean(bev);
bb = mean(bmv);
disp('Best mu and e from taking the average of the best mu and e of each trial')
bestEp2 = sample(floor(bb)) + (sample(floor(bb)+1)-sample(floor(bb)))/(1)*(bb-floor(bb))
bestMu2 = sample(floor(aa)) + (sample(floor(aa)+1)-sample(floor(aa)))/(1)*(aa-floor(aa)) 

%Analysis of Pre-impact cases
vm = v(bmv == 1);
vm2 = v(bmv ~= 1);
l = length(vm);
l2 = length(vm2);
Matrix = zeros(l,7);
Matrix2 = zeros(l2,7);

%Creating dummy variables - these store the coordinates of the places where 
%the radius is close to 0.025 or 0.035
a = [];
a2 = [];
b = [];
b2 = [];
%This information is processed at the very bottom of this code

disp('Statistics for pre-impact data for minimum mu')
for kk = 1:l  
  Matrix(kk,1:6) = bounce_array(vm(kk)).states(1:6);
  d = (bounce_array(vm(kk)).d);   %tangential
  n = (bounce_array(vm(kk)).n);   %normal    
  J = [d;n];
  r1 = sqrt(J(2,3)^2+J(1,3)^2);
  Matrix(kk,7) = r1;
  if r1 < 0.026
      a = [a, vm(kk)];
  elseif r1 > 0.034
      b = [b, vm(kk)];
  end
end
avg = mean(Matrix)
standard = std(Matrix)
Cases = l

disp('Statistics for pre-impact data for all other mu')
for jj = 1:l2
    Matrix2(jj,1:6) = bounce_array(vm2(jj)).states(1:6); 
    d2 = (bounce_array(vm2(jj)).d);   %tangential
    n2 = (bounce_array(vm2(jj)).n);   %normal    
    J2 = [d2;n2];
    r2 = sqrt(J2(2,3)^2+J2(1,3)^2);
    Matrix2(jj,7) = r2;
    if r2 < 0.026
      a2 = [a2, vm2(jj)];
  elseif r2 > 0.034
      b2 = [b2, vm2(jj)];
  end
end
avg2 = mean(Matrix2)
standard2 = std(Matrix2)
Cases = l2

disp('Columns:    x    y    theta    x_dot    y_dot    theta_dot    radius')

%mod(Matrix(:,3),pi/2)*180/pi
%mod(Matrix2(:,3),pi/2)*180/pi
%Creating the plots
tiledlayout(1,2) % Requires R2019b or later
nexttile
histogram(mod(Matrix(:,3),pi)*180/pi,36)
title('Angles for pre-impact data for minimum mu')
nexttile
histogram(mod(Matrix2(:,3),pi)*180/pi,36)
title('Angles for pre-impact data for all other mu')

figure
tiledlayout(1,2) % Requires R2019b or later
nexttile
histogram(Matrix(:,7),40)
title('Radii for pre-impact data for minimum mu')
nexttile
histogram(Matrix2(:,7),40)
title('Radii for pre-impact data for all other mu')

%The following matrices contain some very useful information - they store
%the impact data of those impacts whose radii are very close to matching
%that of one of the major axis
ll1 = length(a);
M1 = zeros(ll1,6);
for hh1 = 1:ll1
    M1(hh1,1:6) = bounce_array(a(hh1)).states(1:6); 
    vect1 = M1(:,3)*180/pi;
end

ll2 = length(b);
M2 = zeros(ll2,6);
for hh2 = 1:ll2
    M2(hh2,1:6) = bounce_array(b(hh2)).states(1:6); 
    vect2 = M2(:,3)*180/pi;
end

ll3 = length(a2);
M3 = zeros(ll3,6);
for hh3 = 1:ll3
    M3(hh3,1:6) = bounce_array(a2(hh3)).states(1:6); 
    vect3 = M3(:,3)*180/pi;
end

ll4 = length(b2);
M4 = zeros(ll4,6);
for hh4 = 1:ll4
    M4(hh4,1:6) = bounce_array(b2(hh4)).states(1:6);
    vect4 = M4(:,3)*180/pi;
end

%Creating a figure where we can display all of this
tiledlayout(2,2) % Requires R2019b or later
nexttile
histogram(vect1,48)
title('Angles for pre-impact data for minimum mu when the radius < 0.026')
nexttile
histogram(vect2,48)
title('Angles for pre-impact data for minimum mu when the radius > 0.034')
nexttile
histogram(vect3,48)
title('Angles for pre-impact data for all other mu when the radius < 0.026')
nexttile
histogram(vect4,48)
title('Angles for pre-impact data for all other mu when the radius > 0.034')
