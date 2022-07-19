% Exercise 2.4

clc;
clear;
% E(1/X) = 1/E(X) ?
n = [10 100 1e3 2e3 5e3 1e4 1e5 1e6];
numberOfExperiments = length(n);

% E[1/X]
m1 = zeros(numberOfExperiments, 1);
% 1/E[X]
m2 = zeros(numberOfExperiments, 1);

% interval
a = 1;
b = 2;

for i = 1:numberOfExperiments
    X = a + (b-a).*rand(n(i), 1);
    m1(i) = mean(1./X);
    m2(i) = 1/mean(X);
end

semilogx(n, m1, "b-*");
hold on;
semilogx(n, m2, "r- .");
xlabel("Number of samples");
legend("E[1/X]", "1/E[X]");

% For a = 1 and b = 2:
% We see that E[1/X] != 1/E[X]
% E[X] = 1.5, 1/E[X] -> 1/1.5 = 0.667 for big n

% For a = 0 and b = 1:
% E[X] = 0.5, 1/E[X] -> 1/0.5 = 2 for big n

% For a = -1 and b = 1:
% E[X] = 0, 1/E[X] -> 1/0 -> oo for big n