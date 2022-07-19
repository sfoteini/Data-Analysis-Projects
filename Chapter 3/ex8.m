% Exercise 3.8

clc;
clear;
close all;

M = 100;        % number of samples
n = 10;         % number of observations in a sample
meanX = 0;
varX = 1;
data = normrnd(meanX, sqrt(varX), n, M);
mu = mean(data);
std = std(data);

% a
% i) 95% confidence interval of std
alpha = 0.05;
ci_std = zeros(M, 2);
for i = 1:M
    [h,p,ci] = vartest(data(:, i), varX, 'Alpha', alpha);
    ci_std(i, :) = sqrt(ci);
end

% ii) 95% confidence interval of std using bootstrap samples
B = 1000;       % number of bootstrap samples
bootstrap_std = bootstrp(B, @std, data);
% Calculate low and upper limit for confidence interval
ci_llimit = floor((B+1)*alpha/2);
ci_ulimit = B+1-ci_llimit;
% Sort the bootstrap_mean array
bootstrap_std_sorted = sort(bootstrap_std);
% Get the confidence interval
bci_std = zeros(M, 2);
for i = 1:M
    bci_std(i,1) = bootstrap_std_sorted(ci_llimit, i);
    bci_std(i,2) = bootstrap_std_sorted(ci_ulimit, i);
end
% Histogram for upper ci limit
figure();
histogram(bci_std(:,1));
hold on;
histogram(ci_std(:,1));
hold off;
title('Histogram of upper ci limit');
legend('Bootstrap ci', 'Parametric ci');
% Histogram for low ci limit
figure();
histogram(bci_std(:,2));
hold on;
histogram(ci_std(:,2));
hold off;
title('Histogram of low ci limit');
legend('Bootstrap ci', 'Parametric ci');

% c : Y = X^2
Ydata = data .* data;
% i) 95% confidence interval of mean
ci_Ystd = zeros(M, 2);
for i = 1:M
    [h,p,ci] = vartest(Ydata(:, i), 1, 'Alpha', alpha);
    ci_Ystd(i, :) = ci;
end
% ii) 95% confidence interval of mean using bootstrap samples
bootstrap_Ystd = bootstrp(B, @std, Ydata);
% Sort the bootstrap_mean array
bootstrap_Ystd_sorted = sort(bootstrap_Ystd);
% Get the confidence interval
bci_Ystd = zeros(M, 2);
for i = 1:M
    bci_Ystd(i,1) = bootstrap_Ystd_sorted(ci_llimit, i);
    bci_Ystd(i,2) = bootstrap_Ystd_sorted(ci_ulimit, i);
end
% Histogram for upper ci limit
figure();
histogram(bci_Ystd(:,1));
hold on;
histogram(ci_Ystd(:,1));
hold off;
title('Histogram of upper ci limit');
legend('Bootstrap ci', 'Parametric ci');
% Histogram for low ci limit
figure();
histogram(bci_Ystd(:,2));
hold on;
histogram(ci_Ystd(:,2));
hold off;
title('Histogram of low ci limit');
legend('Bootstrap ci', 'Parametric ci');