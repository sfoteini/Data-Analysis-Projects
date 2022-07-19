% Exercise 3.5

clc;
clear;
close all;
% Old Faithful eruption data
data = importdata("eruption.dat.txt");

% a : Confidence interval of standard deviation
fprintf("Variance Test\n");
var_waitingtime = 100;
var_duration = 1;
alpha = 0.05;

h_var = zeros(3,1);
p_var = zeros(3,1);
ci_var = zeros(3,2);
% 1st column: waiting time 1989
fprintf("Waiting Time 1989\n");
[h_var(1), p_var(1), ci_var(1,:)] = print_vartest(data(:,1), ...
    var_waitingtime, alpha);
% 2nd column: duration of eruptions 1989
fprintf("\nDuration of eruptions 1989\n");
[h_var(2), p_var(2), ci_var(2,:)] = print_vartest(data(:,2), ...
    var_duration, alpha);
% 3rd column: waiting time 2006
fprintf("\nWaiting Time 2006\n");
[h_var(3), p_var(3), ci_var(3,:)] = print_vartest(data(:,3), ...
    var_waitingtime, alpha);

% b : Confidence interval of mean
fprintf("\nMean value Test\n");
mu_waitingtime = 75;
mu_duration = 2.5;
h_mu = zeros(3,1);
p_mu = zeros(3,1);
ci_mu = zeros(3,2);
% 1st column: waiting time 1989
fprintf("Waiting Time 1989\n");
[h_mu(1), p_mu(1), ci_mu(1,:)] = print_ttest(data(:,1), mu_waitingtime, ...
    alpha);
% 2nd column: duration of eruptions 1989
fprintf("\nDuration of eruptions 1989\n");
[h_mu(2), p_mu(2), ci_mu(2,:)] = print_ttest(data(:,2), mu_duration, ...
    alpha);
% 3rd column: waiting time 2006
fprintf("\nWaiting Time 2006\n");
[h_mu(3), p_mu(3), ci_mu(3,:)] = print_ttest(data(:,3), mu_waitingtime, ...
    alpha);

% c : Goodness-of-fit Test for normal distribution
fprintf("\nGoodness-of-fit Test\n");
h_chi2 = zeros(3,1);
p_chi2 = zeros(3,1);
% 1st column: waiting time 1989
fprintf("Waiting Time 1989\n");
[h_chi2(1), p_chi2(1)] = print_chi2gof(data(:,1), alpha);
figure();
histogram(data(:,1), 'Normalization', 'pdf');
title("Waiting Time (1989)");
% 2nd column: duration of eruptions 1989
fprintf("\nDuration of eruptions 1989\n");
[h_chi2(2), p_chi2(2)] = print_chi2gof(data(:,2), alpha);
figure();
histogram(data(:,2), 'Normalization', 'pdf');
title("Duration of eruptions (1989)");
% 3rd column: waiting time 2006
fprintf("\nWaiting Time 2006\n");
[h_chi2(3), p_chi2(3)] = print_chi2gof(data(:,3), alpha);
figure();
histogram(data(:,3), 'Normalization', 'pdf');
title("Waiting Time (2006)");

% d : Check the claim "With an error of 10 minutes, Old Faithful will erupt 
% 65 minutes after an eruption lasting less than 2.5 minutes or 91 minutes
% after an eruption lasting more than 2.5 minutes."
waitingtimes = data(:, 1);
waitingtime_less = waitingtimes(data(:,2)<2.5);
waitingtime_more = waitingtimes(data(:,2)>=2.5);
mu_less = 65;
mu_more = 91;
std = 10;

fprintf("\nConfidence Interval and null Hypothesis test for mean value\n");
fprintf("Waiting time < 2.5 minutes\n");
[h_mu_l, p_mu_l, ci_mu_l] = print_ttest(waitingtime_less, mu_less, alpha);
fprintf("Waiting time >= 2.5 minutes\n");
[h_mu_m, p_mu_m, ci_mu_m] = print_ttest(waitingtime_more, mu_more, alpha);

fprintf("\nConfidence Interval and null Hypothesis test for std\n");
fprintf("Waiting time < 2.5 minutes\n");
[h_var_l, p_var_l, ci_var_l] = print_vartest(waitingtime_less, std^2, alpha);
fprintf("Waiting time >= 2.5 minutes\n");
[h_var_m, p_var_m, ci_var_m] = print_vartest(waitingtime_more, std^2, alpha);

function [h,p,ci] = print_vartest(data,var,alpha)
    [h,p,ci] = vartest(data, var, 'Alpha', alpha);
    fprintf("%.2f %% Confidence Interval of std: [%.4f, %.4f]\n", ...
        (1-alpha)*100, sqrt(ci(1)), sqrt(ci(2)));
    if h == 1
        fprintf("Rejection of the null Hypothesis H0: σ = %.2f at " + ...
            "%.2f significance level.\n", sqrt(var), alpha);
    else
        fprintf("Acceptance of the null Hypothesis H0: σ = %.2f at " + ...
            "%.2f significance level.\n", sqrt(var), alpha);
    end
end

function [h,p,ci] = print_ttest(data,mu,alpha)
    [h,p,ci] = ttest(data, mu, 'Alpha', alpha);
    fprintf("%.2f %% Confidence Interval of mean value: [%.4f, %.4f]\n", ...
        (1-alpha)*100, ci(1), ci(2));
    if h == 1
        fprintf("Rejection of the null Hypothesis H0: μ = %.2f at " + ...
            "%.2f significance level.\n", mu, alpha);
    else
        fprintf("Acceptance of the null Hypothesis H0: μ = %.2f at " + ...
            "%.2f significance level.\n", mu, alpha);
    end
end

function [h,p] = print_chi2gof(data,alpha)
    [h,p] = chi2gof(data, 'Alpha', alpha);
    if h == 1
        fprintf("Rejection of the null Hypothesis H0: data comes from " + ...
            "a population with a normal distribution at %.2f " + ...
            "significance level.\n", alpha);
    else
        fprintf("Acceptance of the null Hypothesis H0: data comes from " + ...
            "a population with a normal distribution at %.2f " + ...
            "significance level.\n", alpha);
    end
    fprintf("The p-value of the goodness-of-fit test is %.4f.\n", p);
end