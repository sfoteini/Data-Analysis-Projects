% Exercise 3.9

clc;
clear;
close all;

% Data from the X variable
M = 100;        % number of samples
n = 10;         % number of observations in a sample
meanX = 0;
varX = 1;
dataX = normrnd(meanX, sqrt(varX), n, M);
muX = mean(dataX);
stdX = std(dataX);
% Data from the Y variable
M = 100;        % number of samples
m = 12;         % number of observations in a sample
meanY = 0;
% meanY = 0.2;
varY = 1;
dataY = normrnd(meanY, sqrt(varY), m, M);
muY = mean(dataY);
stdY = std(dataY);

% a: Two-sample t-test for equal means
alpha = 0.05;
% i: Parametric 95% Confidence Interval
[h1,p1,ci1] = ttest2(dataX, dataY, 'Alpha', alpha);
% Print confidence interval, null hypothesis test, p-value
%{
for i = 1:M
    fprintf('%.2f %% Confidence Interval of the differnce of 2 means: [%.4f, %.4f].\n', (1-alpha)*100, ci1(1,i), ci1(2,i));
    if h1(i) == 1
        fprintf('Rejection of the null Hypothesis H0: mX = mY at %.2f significance level.\n', alpha);
    else
        fprintf('Acceptance of the null Hypothesis H0: mX = mY at %.2f significance level.\n', alpha);
    end
    fprintf('The p-value is: %.4f.\n', p1(i))
end
%}

figure(1);
subplot(2,1,1);
histogram(ci1(1,:));
title('Histogram of lower limit of conficence interval');
subplot(2,1,2);
histogram(ci1(2,:));
title('Histogram of upper limit of conficence interval');

% ii: Bootstrap 95% Confidence Interval
B = 1000;       % number of bootstrap samples
bootstrap_meanX = bootstrp(B, @mean, dataX);
bootstrap_meanY = bootstrp(B, @mean, dataY);
b_mean = bootstrap_meanX - bootstrap_meanY;
% Calculate low and upper limit for confidence interval
ci_llimit = floor((B+1)*alpha/2);
ci_ulimit = B+1-ci_llimit;
% Sort the bootstrap_mean array
b_mean_sorted = sort(b_mean);
% Get the confidence interval
bci1 = zeros(2, M);
for i = 1:M
    bci1(1, i) = b_mean_sorted(ci_llimit, i);
    bci1(2, i) = b_mean_sorted(ci_ulimit, i);
end

figure(2);
subplot(2,1,1);
histogram(bci1(1,:));
title('Histogram of lower limit of Bootstrap conficence interval');
subplot(2,1,2);
histogram(bci1(2,:));
title('Histogram of upper limit of Bootstrap conficence interval');

% b: Transform Z = X^2
dataX2 = dataX .* dataX;
dataY2 = dataY .* dataY;
% Two-sample t-test for equal means
% i: Parametric 95% Confidence Interval
[~,~,ci2] = ttest2(dataX2, dataY2, 'Alpha', alpha);
figure(3);
subplot(2,1,1);
histogram(ci2(1,:));
title('Histogram of lower limit of conficence interval');
subplot(2,1,2);
histogram(ci2(2,:));
title('Histogram of upper limit of conficence interval');

% ii: Bootstrap 95% Confidence Interval
bootstrap_meanX2 = bootstrp(B, @mean, dataX2);
bootstrap_meanY2 = bootstrp(B, @mean, dataY2);
b_mean2 = bootstrap_meanX2 - bootstrap_meanY2;
% Sort the bootstrap_mean array
b_mean2_sorted = sort(b_mean2);
% Get the confidence interval
bci2 = zeros(2, M);
for i = 1:M
    bci2(1, i) = b_mean2_sorted(ci_llimit, i);
    bci2(2, i) = b_mean2_sorted(ci_ulimit, i);
end

figure(4);
subplot(2,1,1);
histogram(bci2(1,:));
title('Histogram of lower limit of Bootstrap conficence interval');
subplot(2,1,2);
histogram(bci2(2,:));
title('Histogram of upper limit of Bootstrap conficence interval');