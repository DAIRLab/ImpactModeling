function [OptMu, OptEps] = OptimumVarsStandard

ran = randi([1 500],1,80);
MuVec = [];
EpsVec = zeros(1,80);
MuStickVec = [];
stickcount = 0;
for z = 1:80
    n = ran(z);
    [stick,Mu,Ep] = Error(n); %we can eedit output of this code
    if stick == 1
        MuStickVec(end+1) = Mu;
        stickcount = stickcount+1;
        EpsVec(z) = Ep;
    else
        MuVec = [Mu,MuVec];
        EpsVec(z) = Ep;
    end
end

maxstickMu = max(MuStickVec);
MuStick_w = maxstickMu*stickcount; %weigh the value of the max mu of the sticking trials
[i,j] = size(MuVec);
OptMu = (sum(MuVec)+MuStick_w)/(j+stickcount)

[i,j] = size(EpsVec);
OptEps = sum(EpsVec)/j

%Computing standard deviations:
%First, we need to find the vector that contains all the mus:
MuVec2 = zeros(1,stickcount);
MuVec2(1,1:end) = maxstickMu;
MuVec3 = [MuVec,MuVec2];

%Creating the plots
tiledlayout(2,2) % Requires R2019b or later
nexttile
pd = fitdist(EpsVec','Normal');
step = 3/(length(EpsVec)-1);
x = -1:step:2
y = pdf(pd,x);
plot(x,y,'LineWidth',2)
standardEp = pd.sigma;
yl = ylim;
%Making a line at the mean
line([OptEps, OptEps], yl, 'Color', 'blue', 'LineWidth', 3,'LineStyle','--');
%Making a line at 1 standard away
line([OptEps-standardEp,OptEps-standardEp], yl, 'Color', 'cyan', 'LineWidth', 3,'LineStyle','--');
line([0,0], yl, 'Color', 'black', 'LineWidth', 3,'LineStyle','-');
line([OptEps+standardEp,OptEps+standardEp], yl, 'Color', 'cyan', 'LineWidth', 3,'LineStyle','--');
%Making a line at 0 and 1 where our data set ends
line([1,1], yl, 'Color', 'black', 'LineWidth', 3,'LineStyle','-');
ylim([0,1.2])
set(gca,'ytick',[])
title({'Standard Distribution for 80 Random Epsilon', 'Mean = ' num2str(OptEps), 'Standard Deviation = ' num2str(standardEp)})
xlabel('Range of Epsilon')
legend('Standard Curve','Mean','One standard', 'Boundaries');

nexttile
pd2 = fitdist(MuVec3','Normal');
step2 = 3/(length(MuVec3)-1);
x2 = -1:step:2
y2 = pdf(pd2,x2);
plot(x2,y2,'LineWidth',2)
standardMu = pd2.sigma;
av = pd2.mu;
yl2 = ylim;
%Making a line at the mean
line([av, av], yl2, 'Color', 'blue', 'LineWidth', 3,'LineStyle','--');
%Making a line at 1 standard away
line([av-standardMu,av-standardMu], yl2, 'Color', 'cyan', 'LineWidth', 3,'LineStyle','--');
line([1,1], yl2, 'Color', 'black', 'LineWidth', 3,'LineStyle','-');
line([av+standardMu,av+standardMu], yl2, 'Color', 'cyan', 'LineWidth', 3,'LineStyle','--');
%Making a line at 0 and 1 where our data set ends
line([0,0], yl2, 'Color', 'black', 'LineWidth', 3,'LineStyle','-');
ylim([0,1.2])
set(gca,'ytick',[])
title({'Standard Distribution for 80 Random Mu', 'Mean = ' num2str(av), 'Standard Deviation = ' num2str(standardMu)})
xlabel('Range of Mu')
legend('Standard Curve','Mean','One standard', 'Boundaries');

nexttile
histfit(EpsVec,10,'normal')
yl3 = ylim;
%Making a line at the mean
line([OptEps, OptEps], yl3, 'Color', 'blue', 'LineWidth', 3,'LineStyle','--');
%Making a line at 1 standard away
line([OptEps-standardEp,OptEps-standardEp], yl3, 'Color', 'cyan', 'LineWidth', 3,'LineStyle','--');
line([1,1], yl3, 'Color', 'black', 'LineWidth', 3,'LineStyle','-');
line([OptEps+standardEp,OptEps+standardEp], yl3, 'Color', 'cyan', 'LineWidth', 3,'LineStyle','--');
%Making a line at 0 and 1 where our data set ends
line([0,0], yl3, 'Color', 'black', 'LineWidth', 3,'LineStyle','-');
set(gca,'ytick',[])
title('Histogram for 80 Random Epsilon')
xlabel('Range of Epsilon')
xlim([-1,2])
legend('Standard Curve','Mean','One standard', 'Boundaries');

nexttile
histfit(MuVec3,10,'normal')
yl4 = ylim;
%Making a line at the mean
line([av, av], yl4, 'Color', 'blue', 'LineWidth', 3,'LineStyle','--');
%Making a line at 1 standard away
line([av-standardMu,av-standardMu], yl4, 'Color', 'cyan', 'LineWidth', 3,'LineStyle','--');
line([1,1], yl4, 'Color', 'black', 'LineWidth', 3,'LineStyle','-');
line([av+standardMu,av+standardMu], yl4, 'Color', 'cyan', 'LineWidth', 3,'LineStyle','--');
%Making a line at 0 and 1 where our data set ends
line([0,0], yl4, 'Color', 'black', 'LineWidth', 3,'LineStyle','-');
set(gca,'ytick',[])
title('Standard Distribution for 80 Random Mu')
xlabel('Range of Mu')
xlim([-1,2])
legend('Standard Curve','Mean','One standard', 'Boundaries');

set(gcf, 'Position',  [100, 100, 1000, 800])

end