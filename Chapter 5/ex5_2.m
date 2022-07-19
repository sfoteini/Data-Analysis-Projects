% Exercise 5.2

clc;
clear;
close all;

alpha = 0.05;
M = 1000;       % number of samples
n = 20;         % number of observations per sample
L = 1000;       % number of random samples

muX = 0;
muY = 0;
stdX = 1;
stdY = 1;
% correlation coeficient
r = 0;
% r = 0.5;
sigma = [stdX^2 r; r stdY^2];
mu = [muX muY];
X = zeros(n, M);
Y = zeros(n, M);
R = zeros(L+1, M);
for i = 1:M
    tmp_data = mvnrnd(mu, sigma, n);
    X(:,i) = tmp_data(:,1);
    Y(:,i) = tmp_data(:,2);
    % X = X.^2;
    % Y = Y.^2;
    tmp_R = corrcoef(X(:,i), Y(:,i));
    R(1,i) = tmp_R(1,2);

    % Hypothesis test H0: r = 0 using non-parametric test
    for j = 1:L
        tmp_R = corrcoef(X(:,i), Y(randperm(n),i));
        R(j+1,i) = tmp_R(1,2);
    end
end

t = R .* sqrt((n-2) ./ (1-R.^2));
t_sorted = sort(t(2:L+1,:),1);
% low and upper limit
lowlim = round((alpha/2)*L);
upplim = round((1-alpha/2)*L);
tlow = t_sorted(lowlim, :);
tupp = t_sorted(upplim, :);
rej = sum(t(1,:) - tlow < 0 | t(1,:) - tupp > 0);
prej = rej/M;
fprintf("Percentage of rejections of the null hypothesis %.4f%%.\n", prej*100);

% Hypothesis test H0: r = 0 using t-statistic
tc = tinv(1-alpha/2, n-2);
accept = sum(abs(t(1,:))<tc);
fprintf("Percentage of rejections of the null hypothesis %.4f%%.\n", (M-accept)*100/M);