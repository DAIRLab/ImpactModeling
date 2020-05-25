%Randomizer
%What this script does is that on the inner for loop, it selects 80 random
%numbers ranging from 1 to 500 (that is, our test cases), and then it runs
%the Error function (the CompiledCode), delivering for each one of them the
%minimum Mu and Ep - after averaging them, it completes this whole process
%4 more times (the outer for loop), and again averages the result from each
%iteration to deliver the final Mu and Ep

fMu = zeros(1,5);
fEp = zeros(1,5);

for jj = 1:5

%creates a vector with 80 random numbers ranging from 1 to 80
v = randi([1 500],1,80);

Mu = zeros(1,80);
Ep = zeros(1,80);

for k = 1:80
    [bMu, bEp] = Error(v(k));
    Mu(k) = bMu;
    Ep(k) = bEp;
end

finalMu = mean(Mu);
finalEp = mean(Ep);

fMu(jj) = finalMu;
fEp(jj) = finalEp;

end

M = mean(fMu)
E = mean(fEp)

