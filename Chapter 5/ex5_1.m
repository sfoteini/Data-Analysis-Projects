% Exercise 5.1

clc;
clear;
close all;

alpha = 0.05;
M = 1000;       % number of samples
% number of observations per sample
prompt = 'Select 1 for n = 20, 2 for n = 200\n';
choice = input(prompt);
switch choice
    case 1
        n = 20;
    case 2
        % c : Number of observations per sample is 200
        n = 200;
    otherwise
        n = 20;
end
muX = 0;
muY = 0;
stdX = 1;
stdY = 1;
% correlation coeficient
prompt = 'Select 1 for r = 0, 2 for r = 0.5\n';
choice = input(prompt);
switch choice
    case 1
        r = 0;
    case 2
        r = 0.5;
    otherwise
        r = 0;
end
sigma = [stdX^2 r; r stdY^2];
mu = [muX muY];
X = zeros(n, M);
Y = zeros(n, M);

for i = 1:M
    tmp_data = mvnrnd(mu, sigma, n);
    X(:,i) = tmp_data(:,1);
    Y(:,i) = tmp_data(:,2);
end

% d : Z = X^2
prompt = 'Select 4 for question d\n';
choice = input(prompt);
if choice == 4
    X = X.^2;
    Y = Y.^2;
end

% a : 95% confidence interval of correlation coefficient 
% using Fisher transform
R = zeros(M, 1);
RL = zeros(M, 1);
RU = zeros(M, 1);
for i = 1:M
    [tmp_R, ~, tmp_RL, tmp_RU] = corrcoef(X(:,i), Y(:,i), 'Alpha', alpha);
    R(i) = tmp_R(1,2);
    RL(i) = tmp_RL(1,2);
    RU(i) = tmp_RU(1,2);
end
figure();
subplot(2,1,1);
histogram(RL);
title('Histogram of lower limit');
subplot(2,1,2);
histogram(RU);
title('Histogram of upper limit');

% Coeficient interval using Fisher transform
for i = 1:M
    tmp_R = corrcoef(X(:,i), Y(:,i));
    R(i) = tmp_R(1,2);
end
Z = 0.5*log((1+R)./(1-R));
z = norminv(1-alpha/2);
zstd = sqrt(1/(n-3));
ZL = Z - z*zstd;
ZU = Z + z*zstd;
RL = (exp(2*ZL)-1) ./ (exp(2*ZL)+1);
RU = (exp(2*ZU)-1) ./ (exp(2*ZU)+1);
figure();
subplot(2,1,1);
histogram(RL);
title('Histogram of lower limit using Fisher transform');
subplot(2,1,2);
histogram(RU);
title('Histogram of upper limit using Fisher transform');

percentage = 0;
for i = 1:M
    if RL(i) < r & r < RU(i)
        percentage = percentage + 1;
    end
end
fprintf("Το διάστημα εμπιστοσύνης περιέχει την τιμή r = %.2f σε ποσοστό %.4f%%.\n", r, percentage*100/M);

% b : Hypothesis test H0: r = 0 using t-statistic
t = R .* sqrt((n-2) ./ (1-R.^2));
tc = tinv(1-alpha/2, n-2);
accept = sum(abs(t)<tc);
fprintf("Η μηδενική υπόθεση Η0: ρ = 0 απορρίπτεται σε ποσοστό %.4f%%.\n", (M-accept)*100/M);