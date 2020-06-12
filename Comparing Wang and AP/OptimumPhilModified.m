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

%create contour plot
%colormap(flipud(gray))  %match color from Nima paper
%contourf(sample(1:end), sample(1:end), averageError(1:end, 1:end)', 30);
%colorbar('Direction','reverse'); 
%xlabel("Epsilon")
%ylabel("Mu")

%Finding average mu and ep by taking the average of the best of each trial
aa = mean(bev);
bb = mean(bmv);
disp('Best mu and e from taking the average of the best mu and e of each trial')
bestEp2 = sample(floor(bb)) + (sample(floor(bb)+1)-sample(floor(bb)))/(1)*(bb-floor(bb))
bestMu2 = sample(floor(aa)) + (sample(floor(aa)+1)-sample(floor(aa)))/(1)*(aa-floor(aa)) 

%Analysis of Pre-impact cases
vm = find(v(bmv == 1));
vm2 = find(v(bmv ~= 1));
l = length(vm);
l2 = length(vm2);
Matrix = zeros(l,7);
Matrix2 = zeros(l2,7);

disp('Statistics for pre-impact data for minimum mu')
for kk = 1:l  
  Matrix(kk,1:6) = bounce_array(vm(kk)).states(1:6);
  d = (bounce_array(vm(kk)).d);   %tangential
  n = (bounce_array(vm(kk)).n);   %normal    
  J = [d;n];
  Matrix(kk,7) = sqrt(J(2,3)^2+J(1,3)^2);
end
avg = mean(Matrix)
standard = std(Matrix)
Cases = l

disp('Statistics for pre-impact data for all other mu')
for jj = 1:l2
    Matrix2(jj,1:6) = bounce_array(vm2(jj)).states(1:6); 
    d = (bounce_array(vm2(jj)).d);   %tangential
    n = (bounce_array(vm2(jj)).n);   %normal    
    J = [d;n];
    Matrix(jj,7) = sqrt(J(2,3)^2+J(1,3)^2);
end
avg2 = mean(Matrix2)
standard2 = std(Matrix2)
Cases = l2

disp('Columns:    x    y    theta    x_dot    y_dot    theta_dot')

%mod(Matrix(:,3),pi/2)*180/pi
%mod(Matrix2(:,3),pi/2)*180/pi
%Creating the plots
tiledlayout(1,2) % Requires R2019b or later
nexttile
histogram(mod(Matrix(:,3),pi)*180/pi,36)
title('Statistics for pre-impact data for minimum mu')
nexttile
histogram(mod(Matrix2(:,3),pi)*180/pi,36)
title('Statistics for pre-impact data for all other mu')
