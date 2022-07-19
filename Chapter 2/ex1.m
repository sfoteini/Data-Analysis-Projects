% Exercise 2.1

clc;
clear;
% Coin flipping experiment
nOfRepetitions = [10 100 1e3 2e3 5e3 1e4 1e5 1e6];
numberOfExperiments = length(nOfRepetitions);
tailsFrequency = zeros(6,1);

% experiments
for i=1 : numberOfExperiments
    flips = unidrnd(2,nOfRepetitions(i),1);
    nOfTails = length(flips(flips==1));
    tailsFrequency(i) = (nOfTails /nOfRepetitions(i)) * 100;
    fprintf('%d coin flips: %f %%.\n', nOfRepetitions(i),tailsFrequency(i));
end
figure();
plot(log10(nOfRepetitions), tailsFrequency);
xlabel('Number of repetitions (log10)');
ylabel('Probaility of tails');
title('Probability of a binary outcome');