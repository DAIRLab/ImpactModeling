function M = explore_angles(bounce_array,trial)

%trial number to draw
trial = 700;

%set up pre-impact ellipse
ang = bounce_array(trial).states(3);
x0 = bounce_array(trial).states(1);
y0 = bounce_array(trial).states(2);

[xvec,yvec] = ellipse_vector(ang,x0,y0);
pos = find(yvec == min(yvec));
y1 = min(yvec);
x1 = xvec(pos);

%variation points
range = ang*0.1;
sample_ang = linspace(ang-range, ang+range, 51);
x1_vec = zeros(1,51);
y1_vec = zeros(1,51);

for i = 1:51
    ang2 = sample_ang(i);
    [xvec3,yvec3] = ellipse_vector(ang2,x0,y0);
    pos3 = find(yvec3 == min(yvec3));
    x1_vec(i) = min(xvec3(pos3));
    y1_vec(i) = min(yvec3);
end

%split between positive and negative y values. We are only interested in
%the positive y values
a = find(y1_vec >= 0);
b = find(y1_vec < 0);

%y1_positive = y1_vec(a);
%x1_positive = x1_vec(a);
y1_positive = y1_vec(a);
x1_positive = x1_vec(a);

%begin a double for loop to test each of the values in the y/x positive
%vectors

% Step 2: set-up the interval
iter = 50; %how many iterations of mu and epsilon we would like
sample = linspace(0.01, 0.99, iter); %set up sample vector
trialError = zeros(iter);
totalError = zeros(iter);

% Step 3: run the simulation
%Get pre and post impact data from current trial
pre = bounce_array(trial).states(4:6);
post = bounce_array(trial).states(10:12);

%get normal and tangenial Jacobians from dataset
d = (bounce_array(trial).d);   %tangential
n = (bounce_array(trial).n);   %normal

y1_positive(end+1) = n(1,3);
x1_positive(end+1) = d(1,3);

J = [d;n];

%establishing the percentages
p1 = 1;
p2 = 1;

M = zeros(3,length(y1_positive));

for l = 1:length(y1_positive)
    
    %modifying J
    d(1,3) = x1_positive(l);
    n(1,3) = y1_positive(l);

    %loop over all mu's
    for i = 1:iter
        %assign epsilon to value corresponding with loop
        epsilon = sample(i);

        %loop over all epsilons
        for j = 1:iter
            %assign mu to value corresponding with loop
            mu = sample(j);

            %Run AP Poisson Model given mu and epsilon
             %v_calc = APPoisson_juniors(M, n, d, pre, mu, epsilon);
             v_calc = Wang_juniors(pre,n,d,mu,epsilon,p1,p2);

             error = sqrt((post(1) - v_calc(1))^2 + (post(2) - v_calc(2))^2)/sqrt(post(1)^2 + post(2)^2); 

             trialError(i, j) = error;            

        end

    end
    
    % Step 4: post-process results
    %Finding average mu and ep by taking the average of the best of each trial
    minError = min(min(trialError));
    [be , bm] = find(trialError == min(min(trialError)));
    M(1,l) = minError;
    M(2,l) = min(sample(be));
    M(3,l) = min(sample(bm));
    if l < length(y1_positive)
        M(4,l) = sample_ang(a(l));
    else
        M(4,l) = ang;
    end
end

end
    

