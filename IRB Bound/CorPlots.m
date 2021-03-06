%% Test Some Ideas

difference = abs(noWidth - yesWidth(1, :));
plot(difference, yesWidth(2, :), '.')
ylabel("Optimal Width [m]")
xlabel("Difference In Omega Dot Accuracy")
xlim([0, 0.05])


%% 
count = 0;
for i = 1:200
    trial = i;
    if sum(bounce_array(trial).flags) < 1
        count = count + 1;
        pre = bounce_array(trial).states(4:6)';
        post = bounce_array(trial).states(10:12)';
        theta = (rem(bounce_array(trial).states(3),pi)*180)/pi;
        out(1, count) = theta;
        out(2, count) = post(3) - pre(3);
    end
end

hold on
plot(out(1,:), out(2,:), '.')
plot(linspace(-180, 180), -80*sind(2*linspace(-180,180)))
ylim([-90, 90]);

xlabel("Pre Impact Theta [deg]")
ylabel("Change in Rotational Velocity (Observed) [rad/s]")

%%
figure()
hold on
plot(out(1,:), checkNoWidth(2,:), '.')
%plot(out(1,:), checkWidth(2,:), '.')
%plot(ones(100,1)*45, linspace(0, 1.5));
xlabel("Pre-Impact \theta [deg]")
ylabel("Error")
title("IRB (Without Width)");

%% 
figure();
plot(checkWidth(3, :), bestWidth(3,:), '.')
xlabel("Net Moment");
ylabel("Optimal Width");

%%
figure()
hold on
plot(1:173, checkWith(:), '.');
plot(1:173, checkNoWidth(:), 'r.');
ylim([-0.04, 0.04]);
xlabel("Trial Number")
ylabel("|\omega_{observed}| - |\omega_{predicted}|")
legend("IRB With Patch Width", "IRB Without Patch Width", "Location", "North")
