% Exercise 3.1

clc;
clear;
% Poisson distribution
% a : Show that the maximum likelihood estimator of lambda is the 
% mean value
% The likelihood function
% L(x1, ..., xn; λ) = Π exp(-λ)*λ^(xj)/(xj!)
% The log-likelihood function
% logL(x1, ..., xn; λ) = -nλ - Σln(xj!) + ln(λ) * Σxj
% dlogL/dλ = 0 -> -n + Σxj / λ = 0 -> λ = Σxj / n

n = [10 1e2 1e3 1e4 1e5 1e6];
lambda = 3;
fprintf('Lambda: %d\n', lambda);
meanX = zeros(length(n),1);
for i = 1:length(n)
    X = poissrnd(lambda, n(i), 1);
    meanX(i) = mean(X);
    fprintf('Number of samples: %d Mean: %.5f\n', n(i), meanX(i));
end

% b
m1 = poissonMean(1, 10, 100);
m2 = poissonMean(2, 100, 1e4);
m3 = poissonMean(3, 100, 1e4);


function m = poissonMean(lambda, M, n)
   X = poissrnd(lambda, n, M);
   meanX = mean(X);
   figure();
   histogram(meanX);
   title(sprintf('Mean values from poisson distribution lambda=%.2f', ...
       lambda));
   m = mean(meanX);
   fprintf('Mean value of Poisson mean values: %.5f\n', m);
   hold on;
   xline(m);
   hold off;
end