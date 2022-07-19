% Exercise 3.6

clc;
clear;
close all;
% Bootstrap estimate of standard error of sample mean
n = 10;
meanX = 0;
varX = 1;
data = normrnd(meanX, sqrt(varX), 1, n);
mu = mean(data);

% a : Mean of bootstrap samples
B = 1000;       % number of bootstrap samples
bootstap_samples = zeros(n, B);
for i = 1:B
    tmp_index = randi(n, 1, n);
    tmp_data = data(tmp_index);
    bootstap_samples(:, i) = tmp_data;
end
bsamples_mean = mean(bootstap_samples);
figure();
histogram(bsamples_mean);
xline(mu, 'Color', 'r', 'LineWidth', 2);
xline(meanX, 'Color', 'g', 'LineWidth', 2, 'LineStyle', '--');
legend('bootstrap mean', 'data mean', 'X mean');
title('Mean value of bootstrap samples');
xlabel('Mean value');
ylabel('Count');
mean_bootstrap = mean(bsamples_mean);
fprintf('Mean value from data m = %.3f\n', mu);
fprintf('Mean value from bootstrap samles m = %.3f\n', mean_bootstrap);

% b : Calculate standard error
serror = std(bsamples_mean);
serror_data = std(data)/sqrt(n);
fprintf('Standard Error from data se = %.3f\n', serror_data);
fprintf('Standard Error from bootstrap samles se = %.3f\n', serror);

% c : Y = exp(X)
Ydata = exp(data);
muY = mean(Ydata);
Ybootstap_samples = zeros(n, B);
for i = 1:B
    tmp_index = randi(n, 1, n);
    tmp_data = Ydata(tmp_index);
    Ybootstap_samples(:, i) = tmp_data;
end
Ybsamples_mean = mean(Ybootstap_samples);
figure();
histogram(Ybsamples_mean);
xline(muY, 'Color', 'r', 'LineWidth', 2);
legend('bootstrap mean', 'data mean');
title('Mean value of bootstrap samples');
xlabel('Mean value');
ylabel('Count');
meanY_bootstrap = mean(Ybsamples_mean);
fprintf('Mean value from data m = %.3f\n', muY);
fprintf('Mean value from bootstrap samles m = %.3f\n', meanY_bootstrap);

Yserror = std(Ybsamples_mean);
Yserror_data = std(Ydata)/sqrt(n);
fprintf('Standard Error from data se = %.3f\n', Yserror_data);
fprintf('Standard Error from bootstrap samles se = %.3f\n', Yserror);