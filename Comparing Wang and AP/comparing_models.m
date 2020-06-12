% Comparing Wang, AP script

% load in ellipse data 
load('ellipse_uniform.mat');

% Set up Constants
a0 = 0.07/2; %semi-major axis
b0 = 0.05/2; %semi-minor axis
m1 = 0.037; 
I1 = m1 * (a0^2 + b0^2) / 4;
M = [m1,0,0;0,m1,0;0,0,I1];
AP_g = 0;
Wang_g = 0;
AP_i = 0;
Wang_i = 0;


% Set up interval
iter = 50; %how many iterations of mu and epsilon we would like
sample = linspace(0.01, 0.99, iter); %set up sample vector
trialError = zeros(iter);
models = zeros(iter);
totalError = zeros(iter);

% Run Simulation
%loop over multiple trials
numTrials = 2000; %Number of Trials
bev = zeros(1,numTrials);
bmv = zeros(1,numTrials);
v = randi([1 numTrials],1,2000);

% for each Trial, call both the Wang and AP models,
% 1) receive from each of them the optimal mu and e for that case,
% as well as the post-impact velocities
% 2) simultaneously compute the associated error from the 
% actual post-impact velocity
% 3) choose the "best" model of the two in each trial, and create
% a growing vector with the best mu and e

for trial = 1:numTrials
    %trial = v(tr);
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

            %Run AP Poisson Model given mu and epsilon, compute error 
             v_calc_ap = APPoisson_juniors(M, n, d, pre, mu, epsilon);
             
             %Pt_ap = (Ma*(J*(post-pre)'));
             %p_hat_ap = (Ma*(J*(v_calc_ap-pre')));
             %error_ap = norm(Pt_ap-p_hat_ap);
             %error_ap = error_ap/norm(Pt_ap);
             error_ap = sqrt((post(1) - v_calc_ap(1))^2 + (post(2) - v_calc_ap(2))^2)/sqrt(post(1)^2 + post(2)^2); 
             %error_ap = 5;
             
            %Run Wang Model given mu and epsilon, compute error
             v_calc_wang = Wang_juniors(pre,n,d,mu,epsilon);
             
             %Pt_wang = (Ma*(J*(post-pre)'));
             %p_hat_wang = (Ma*(J*(v_calc_wang-pre')));            
             %error_wang = norm(Pt_wang-p_hat_wang);
             %error_wang = error_wang/norm(Pt_wang);
             error_wang = sqrt((post(1) - v_calc_wang(1))^2 + (post(2) - v_calc_wang(2))^2)/sqrt(post(1)^2 + post(2)^2); 
            
             if error_wang < error_ap
                 trialError(i, j) = error_wang;  
                 models(i, j) = 1; 
             else
                 trialError(i, j) = error_ap; 
                 models(i, j) = 2; 
             end
            
        end
    end
    
    [bm , be] = find(trialError == min(min(trialError)));
    Wang_g = Wang_g + sum(sum(models == 1));
    AP_g = AP_g + sum(sum(models == 2));
    if models(bm,be) == 1
        Wang_i = Wang_i + 1;
    else
        AP_i = AP_i + 1;
    end
    bev(tr) = min(be);
    bmv(tr) = min(bm);
    
    %add trial error to total error
    totalError = totalError +  trialError;
    
end


%% Post Process Results
%Take the average of totalError
averageError = totalError / numTrials;
minError = min(min(averageError));  %get the minimum error
[a , b] = find(averageError == minError); %find its index
%convert the index to mu and epsilon values
disp('Best mu and e from taking the average of the added matrices')
bestEp = sample(a)
bestMu = sample(b)

%create contour plot
colormap(flipud(gray))  %match color from Nima paper
contourf(sample(1:end), sample(1:end), averageError(1:end, 1:end)', 20);
colorbar('Direction','reverse'); 

xlabel("Epsilon")
ylabel("Mu")

%Computing average mu and e
aa = mean(bev);
bb = mean(bmv);
disp('Best mu and e from taking the average of the best mu and e of each trial')
bestEp2 = sample(floor(bb)) + (sample(floor(bb)+1)-sample(floor(bb)))/(1)*(bb-floor(bb))
bestMu2 = sample(floor(aa)) + (sample(floor(aa)+1)-sample(floor(aa)))/(1)*(aa-floor(aa)) 
disp('Percentage of times Wang and AP are chosen in the error matrix')
Wang_gper = Wang_g/(Wang_g+AP_g)
AP_gper = AP_g/(Wang_g+AP_g)
disp('Percentage of times Wang and AP are chosen for best mu and e')
Wang_iper = Wang_i/(Wang_i+AP_i)
AP_iper = AP_i/(Wang_i+AP_i)
