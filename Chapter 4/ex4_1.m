% Exercise 4.1

clc;
clear;
close all;

% a
h1 = 100;
h2 = [60 54 58 60 56];
alpha = 0.05;
n = length(h2);
% coefficient of restitution
e0 = 0.76;
e = sqrt(h2/h1);
std_e = std(e);
t = tinv(1-alpha/2, n-1);
mu_e = mean(e);
sderr_e = std_e/sqrt(n);
% accuracy mean(e) - mean(e0)
fprintf('Precision of COR: %.4f\n', std_e);
fprintf('Precision limit of COR: e +- t1-a/2, n-1 * s = e +- %.4f\n', t*std_e);
fprintf('Precision of mean COR: %.4f\n', sderr_e);
fprintf('Precision limit of mean COR: mean{e} +- t1-a/2, n-1 * std error = %.4f +- %.4f\n', mu_e, t*sderr_e);

% b
M = 1000;       % number of samples
n = 5;          % number of observations
mu = 58;
sigma = 2;
h2 = normrnd(mu, sigma, n, M);
mu_h2 = mean(h2);
std_h2 = std(h2);
% coefficient of restitution
e = zeros(n, M);
for i = 1:M
    e(:,i) = sqrt(h2(:,i)./h1);
end
mu_e = mean(e);
std_e = std(e);
emean = sqrt(mu/h1);
esigma = 0.5*sqrt(1/h1)*sqrt(1/mu)*sigma;

figure();
histogram(mu_h2);
title('Histogram of mean values of h2');
xline(mu, 'Color', 'r', 'LineWidth', 2);
legend('mean of h2', 'mean');
figure();
histogram(std_h2);
title('Histogram of std of h2');
xline(sigma, 'Color', 'r', 'LineWidth', 2);
legend('std of h2', 'sigma');
figure();
histogram(mu_e);
title('Histogram of mean values of e');
xline(emean, 'Color', 'r', 'LineWidth', 2);
legend('mean of e', 'mean');
figure();
histogram(std_e);
title('Histogram of std of e');
xline(esigma, 'Color', 'r', 'LineWidth', 2);
legend('std of e', 'sigma');

% c
h1 = [80 100 90 120 95];
h2 = [48 60 50 75 56];
mu_h1 = mean(h1);
std_h1 = std(h1);
mu_h2 = mean(h2);
std_h2 = std(h2);
e = sqrt(h2./h1);
mu_e = mean(e);
std_e = std(e);

% h1, h2 independent
esigma = sqrt((-0.5*sqrt(mu_h2)*mu_h1^(-3/2))^2*std_h1^2 + ...
    (0.5*sqrt(1/mu_h1)*sqrt(1/mu_h2))^2*std_h2^2);
fprintf('Uncertainty of COR: %.4f\n', std_e);
fprintf('Propagation of error of COR: %.4f\n', esigma);

% h1, h2 dependent to each other
r =  corrcoef(h1,h2);
r = r(1,2);
h1h2covar = std_h1 * std_h2 * r;
esigma2 = sqrt((-0.5*sqrt(mu_h2)*mu_h1^(-3/2))^2*std_h1^2 + ...
    (0.5*sqrt(1/mu_h1)*sqrt(1/mu_h2))^2*std_h2^2 + ...
    2*(-0.5*sqrt(mu_h2)*mu_h1^(-3/2))*(0.5*sqrt(1/mu_h1)* ...
    sqrt(1/mu_h2))*h1h2covar);
fprintf('Propagation of error of COR: %.4f\n', esigma2);

% Accuracy mean(e) - e0 -> systematic error
% Random error -> precision boundaries mean +- t * s