% Exercise 2.5

clc;
clear;
% Probabilities on normal distribution
% length: X ~ N(4, 0.01)
length_limit = 3.9;
mean = 4;
sigma = 0.1;
% probability: P(X < 3.9) = Φ((3.9-4)/0.1)
pr_destroy = normcdf((length_limit - mean)/sigma);
fprintf("The probability of a rail being destroyed is %f or %.2f%%.\n", pr_destroy, pr_destroy*100);
% probability: P(X < a) = 0.01 -> Φ((a-4)/0.1) = 0.01
p = 0.01;
length = norminv(p,mean,sigma);
fprintf("The lower limit of length is %f.\n", length);
