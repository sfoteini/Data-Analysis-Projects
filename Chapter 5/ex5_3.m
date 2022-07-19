% Exercise 5.3

clc;
clear;
close all;

rain = importdata("rainThes59_97.dat.txt");
temp = importdata("tempThes59_97.dat.txt");
[y, m] = size(rain);        % y : years, number of observations per sample
                            % m : months, number of samples
L = 1000;                   % number of random permutations per sample
alpha = 0.05;
R = zeros(L+1, m);

for i = 1:m
    tmp_R = corrcoef(rain(:,i), temp(:,i));
    R(1,i) = tmp_R(1,2);
    % Hypothesis test H0: r = 0 using non-parametric test
    for j = 1:L
        tmp_R = corrcoef(rain(:,i), temp(randperm(y),i));
        R(j+1,i) = tmp_R(1,2);
    end
end

t = R .* sqrt((y-2) ./ (1-R.^2));
t_sorted = sort(t(2:L+1,:),1);
% low and upper limit
lowlim = round((alpha/2)*L);
upplim = round((1-alpha/2)*L);
tlow = t_sorted(lowlim, :);
tupp = t_sorted(upplim, :);
rej = t(1,:) - tlow < 0 | t(1,:) - tupp > 0;
nrej = sum(rej);
prej = nrej/m;
fprintf("Percentage of rejections of the null hypothesis %.4f%%.\n", prej*100);

% Hypothesis test H0: r = 0 using t-statistic
tc = tinv(1-alpha/2, y-2);
accept = abs(t(1,:))<tc;
naccept = sum(accept);
fprintf("Percentage of rejections of the null hypothesis %.4f%%.\n", (m-naccept)*100/m);

% Print stats per month
for i = 1:m
    fprintf("--------- Month %d ---------\n", i);
    % Non-parametric
    fprintf("Non-parametric: t = %.3f, ci = [%.3f, %.3f] ", t(1,i), tlow(1,i), tupp(1,i));
    if rej(1,i) == 1
        fprintf("Rejection\n");
    else
        fprintf("Acceptance\n");
    end
    % Parametric
    fprintf("Parametric: t = %.3f, tinv = %.3f ", t(1,i), tc);
    if accept(1,i) == 0
        fprintf("Rejection\n");
    else
        fprintf("Acceptance\n");
    end
end