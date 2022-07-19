% Exercise 3.10

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
alpha = [0.025 0.05 0.075];
for a = 1:length(alpha)
    % i: Parametric
    % Null hypothesis H0: mX = 0
    [h1,p1] = ttest(data, 0, 'Alpha', alpha(a));
    % Null hypothesis H0: mX = 0.5
    [h2,p2] = ttest(data, 0.5, 'Alpha', alpha(a));
    
    figure();
    subplot(4,2,1);
    histogram(p1);
    title("Histogram of p-value for H0: m = 0");
    subplot(4,2,2);
    histogram(p2);
    title("Histogram of p-value for H0: m = 0.5");
    subplot(4,2,3);
    histogram(h1,2);
    title("Null Hypothesis H0: m = 0");
    subplot(4,2,4);
    histogram(h2,2);
    title("Null Hypothesis H0: m = 0.5");
    fprintf("PARAMETRIC\n")
    fprintf("Percentage of rejections for H0: m = 0 at %.4f significance " + ...
        "level: %.3f%%.\n", alpha(a), length(h1(h1==1))/M*100);
    fprintf("Percentage of rejections for H0: m = 0.5 at %.4f significance " + ...
        "level: %.3f%%.\n", alpha(a), length(h2(h2==1))/M*100);
    
    % ii: Bootstrap
    m0 = [0 0.5];
    B = 1000;       % Number of bootstrap samples
    bootstrap_mean = zeros(B+1, M);
    bootstrap_p = zeros(M,1);
    fprintf("BOOTSTRAP\n")
    for i = 1:length(m0)
        data_norm = data - mu + m0(i);
        bootstrap_mean(1:B,:) = bootstrp(B, @mean, data_norm);
        bootstrap_mean(B+1,:) = mu;
        bootstrap_mean_sorted = sort(bootstrap_mean);
        rejection = 0;
        for j = 1: M
            r = find(bootstrap_mean_sorted(:,j) == mu(j));      %rank
            % if all the values are identical, select the middle rank
            if length(r) == B+1
                r = round((B+1)/2);
            elseif length(r) >= 2
                % If at least one bootstrap statistic is identical to the 
                % original, pick the rank of one of them at random
                r = r(unidrnd(length(r)));
            end
            if r > 0.5*(B+1)
                bootstrap_p(j) = 2*(1-r/(B+1));
            else
                bootstrap_p(j) = 2*r/(B+1);
            end   
            if r < (B+1)*alpha(a)/2 | r > (B+1)*(1-alpha(a)/2)
                rejection = rejection + 1;
            end
        end
        fprintf("Percentage of rejections for H0: m = %.2f at %.4f " + ...
            "significance level: %.3f%%.\n", m0(i), alpha(a), ...
            rejection/M*100);
        figure();
        histogram(bootstrap_p);
    end
    fprintf("\n");
end

% b: Transform Y = X^2
dataY = data .* data;
muY = mean(data);

for a = 1:length(alpha)
    % i: Parametric
    % Null hypothesis H0: mY = 1
    [h1,p1] = ttest(dataY, 1, 'Alpha', alpha(a));
    % Null hypothesis H0: mY = 2
    [h2,p2] = ttest(dataY, 2, 'Alpha', alpha(a));
    
    figure();
    subplot(4,2,1);
    histogram(p1);
    title("Histogram of p-value for H0: m = 1");
    subplot(4,2,2);
    histogram(p2);
    title("Histogram of p-value for H0: m = 2");
    subplot(4,2,3);
    histogram(h1,2);
    title("Null Hypothesis H0: m = 1");
    subplot(4,2,4);
    histogram(h2,2);
    title("Null Hypothesis H0: m = 2");
    fprintf("PARAMETRIC\n")
    fprintf("Percentage of rejections for H0: m = 1 at %.4f significance level: %.3f%%.\n", alpha(a), length(h1(h1==1))/M*100);
    fprintf("Percentage of rejections for H0: m = 2 at %.4f significance level: %.3f%%.\n", alpha(a), length(h2(h2==1))/M*100);
    
    % ii: Bootstrap
    m0 = [1 2];
    bootstrap_mean = zeros(B+1, M);
    fprintf("BOOTSTRAP\n")
    for i = 1:length(m0)
        data_norm = dataY - muY + m0(i);
        bootstrap_mean(1:B,:) = bootstrp(B, @mean, data_norm);
        bootstrap_mean(B+1,:) = muY;
        bootstrap_mean_sorted = sort(bootstrap_mean);
        rejection = 0;
        for j = 1: M
            r = find(bootstrap_mean_sorted(:,j) == muY(j));      %rank
            if r < (B+1)*alpha(a)/2 | r > (B+1)*(1-alpha(a)/2)
                rejection = rejection + 1;
            end
        end
        fprintf("Percentage of rejections for H0: m = %.2f at %.4f significance level: %.3f%%.\n", m0(i), alpha(a), rejection/M*100);
    end
    fprintf("\n");
end
