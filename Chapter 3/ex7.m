% Exercise 3.7

clc;
clear;
close all;

M = 100;        % number of samples
n = 10;         % number of observations in a sample
meanX = 0;
varX = 1;
data = normrnd(meanX, sqrt(varX), n, M);
mu = mean(data);

% a
% i) 95% confidence interval of mean
alpha = 0.05;
ci_mean = zeros(M, 2);
for i = 1:M
    [h,p,ci] = ttest(data(:, i), meanX, 'Alpha', alpha);
    ci_mean(i, :) = ci';
end

% ii) 95% confidence interval of mean using bootstrap samples
B = 1000;       % number of bootstrap samples
bootstrap_mean = bootstrp(B, @mean, data);
% Calculate low and upper limit for confidence interval
ci_llimit = floor((B+1)*alpha/2);
ci_ulimit = B+1-ci_llimit;
% Sort the bootstrap_mean array
bootstrap_mean_sorted = sort(bootstrap_mean);
% Get the confidence interval
bci_mean = zeros(M, 2);
for i = 1:M
    bci_mean(i,1) = bootstrap_mean_sorted(ci_llimit, i);
    bci_mean(i,2) = bootstrap_mean_sorted(ci_ulimit, i);
end
% Histogram for upper ci limit
figure();
histogram(bci_mean(:,1));
hold on;
histogram(ci_mean(:,1));
hold off;
title('Histogram of upper ci limit');
legend('Bootstrap ci', 'Parametric ci');
% Histogram for low ci limit
figure();
histogram(bci_mean(:,2));
hold on;
histogram(ci_mean(:,2));
hold off;
title('Histogram of low ci limit');
legend('Bootstrap ci', 'Parametric ci');

% c : Y = X^2
Ydata = data .* data;
% i) 95% confidence interval of mean
ci_Ymean = zeros(M, 2);
for i = 1:M
    [h,p,ci] = ttest(Ydata(:, i), 0, 'Alpha', alpha);
    ci_Ymean(i, :) = ci;
end
% ii) 95% confidence interval of mean using bootstrap samples
bootstrap_Ymean = bootstrp(B, @mean, Ydata);
% Sort the bootstrap_mean array
bootstrap_Ymean_sorted = sort(bootstrap_Ymean);
% Get the confidence interval
bci_Ymean = zeros(M, 2);
for i = 1:M
    bci_Ymean(i,1) = bootstrap_Ymean_sorted(ci_llimit, i);
    bci_Ymean(i,2) = bootstrap_Ymean_sorted(ci_ulimit, i);
end
% Histogram for upper ci limit
figure();
histogram(bci_Ymean(:,1));
hold on;
histogram(ci_Ymean(:,1));
hold off;
title('Histogram of upper ci limit');
legend('Bootstrap ci', 'Parametric ci');
% Histogram for low ci limit
figure();
histogram(bci_Ymean(:,2));
hold on;
histogram(ci_Ymean(:,2));
hold off;
title('Histogram of low ci limit');
legend('Bootstrap ci', 'Parametric ci');
