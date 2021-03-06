%% Statiscal Anlysis of Square Data w/ Wang Mason Model

%Phil

%% Examing Mean Square Error 
%Run through all simulations and get best Mus and Epsilons for each trial
numTrials = 500; %total number of trials
cur = 1; %iterator variable

MuVec = [];
EpsVec = [];
MuStickVec = [];


for i = 1:numTrials
        [stick, Mu, Ep] = Error(i); %find optimal parameters for a single trial
        if stick == 1
            MuStickVec(end+1) = Mu;
            %stickcount = stickcount+1;
            EpsVec(i) = Ep;
        else
            MuVec(cur) = Mu;
            EpsVec(i) = Ep;
            cur = cur + 1;
        end
        if Mu == 0.99
            disp("wtf " + i)
        end
end


% create bar graph
barV = zeros(1,20);
for i = 1:length(MuVec)
    Indx = ceil(20*MuVec(i));

    barV(Indx) = barV(Indx) + 1;
end 

bar(0.05:0.05:1,barV);
disp(length(MuVec));
title("Optimal Mu for All NonSticking Trials")
xlabel("Mu")
ylabel("Frequency")



