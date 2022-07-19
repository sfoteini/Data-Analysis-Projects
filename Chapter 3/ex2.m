% Exercise 3.2

clc;
clear;
close all;
% Exponential distribution
% a : Show that the maximum likelihood estimator of 1/lambda is the 
% mean value
% The likelihood function
% L(x1, ..., xn; λ) = Π λ*exp(-λxj)
% The log-likelihood function
% logL(x1, ..., xn; λ) = Σln(λ*exp(-λxj)) = n*ln(λ)-λ*Σxj
% dlogL/dλ = 0 -> n/λ - Σxj = 0 -> 1/λ = Σxj / n

n = [10 1e2 1e3 1e4 1e5 1e6];
lambda = 1/5;
fprintf('Lambda: %d\n', lambda);
meanX = zeros(length(n),1);
for i = 1:length(n)
    X = exprnd(1/lambda, n(i), 1);
    meanX(i) = mean(X);
    fprintf('Number of samples: %d Mean: %.5f\n', n(i), meanX(i));
end

% b
m1 = exponentialMean(1, 10, 100);
m2 = exponentialMean(2, 100, 1e4);
m3 = exponentialMean(3, 100, 1e6);


function m = exponentialMean(lambda, M, n)
   X = exprnd(1/lambda, n, M);
   meanX = mean(X);
   figure();
   histogram(meanX);
   title(sprintf(['Mean values from exponential distribution ' ...
       'lambda=%.2f'],lambda));
   m = mean(meanX);
   fprintf('Mean value of Exponential mean values: %.5f\n', m);
   hold on;
   xline(m);
   hold off;
end