%%
clear;
load('ellipse_uniform.mat');

trials = 2000;
count = 0;
for i = 1:trials
    trial = i;
    if sum(bounce_array(trial).flags) < 1
        count = count + 1;
        [P,w, error] = IRB_NEW_torque(trial);
        errorVec(count) = error;       
    end
end
%%
errorVec(819) = 0;
errorVec(1420) = 0;
errorVec(550) = 0;
errorVec(1419) = 0;

plot(1:count, errorVec, '.')
xlabel("Trail #");
ylabel("l2 Norm Error of Velocity (not normalized)")


    