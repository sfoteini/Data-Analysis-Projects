% Exercise 3.11

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
m = 12;         % number of observations in a sample
meanY = 0;
varY = 1;
dataY = normrnd(meanY, sqrt(varY), m, M);
muY = mean(dataY);
stdY = std(dataY);

B = 1000;       % Number of Bootstrap samples
data = [dataX; dataY];
b_dif_means = zeros(B,100);
bootstap_samplesX = zeros(n,M);
bootstap_samplesY = zeros(m,M);
for i=1:B
    tmp_index = randi(n+m, n+m, 100);
    tmp_data = zeros(n+m, M);
    for j=1:M
        tmp_data(:,j) = data(tmp_index(:,j));
        bootstap_samplesX(:,j) = tmp_data(1:n,j);
        bootstap_samplesY(:,j) = tmp_data(n+1:n+m,j);
    end
    b_dif_means(i,:) = mean(bootstap_samplesX) - mean(bootstap_samplesY);
end
b_dif_means(B+1,:) = muX-muY;
b_dif_means = sort(b_dif_means);

rp_dif_means = zeros(B,100);
rp_samplesX = zeros(n,M);
rp_samplesY = zeros(m,M);
for i=1:B
    [~, tmp_index] = sort(rand(n+m, M));
    tmp_data = zeros(n+m, M);
    for j=1:M
        tmp_data(:,j) = data(tmp_index(:,j));
        rp_samplesX(:,j) = tmp_data(1:n,j);
        rp_samplesY(:,j) = tmp_data(n+1:n+m,j);
    end
    rp_dif_means(i,:) = mean(rp_samplesX) - mean(rp_samplesY);
end
rp_dif_means(B+1,:) = muX-muY;
rp_dif_means = sort(rp_dif_means);

alpha = [0.025 0.05 0.075];
for a = 1:length(alpha)
    fprintf("PARAMETRIC\n")
    h_param = ttest2(dataX, dataY, 'Alpha', alpha(a));
    fprintf("Percentage of rejections for H0: mX = mY at %.4f significance level: %.3f%%.\n", alpha(a), length(h_param(h_param==1))/M*100);
    
    fprintf("BOOTSTRAP\n");
    b_rejections = 0;
    for j=1:M
        r = find(b_dif_means(:,j) == muX(j)-muY(j));      %rank
        % if all the values are identical, select the middle rank
        if length(r) == B+1
            r = round((B+1)/2);
        elseif length(r) >= 2
            % If at least one bootstrap statistic is identical to the 
            % original, pick the rank of one of them at random
            r = r(unidrnd(length(r)));
        end
        if r < (B+1)*alpha(a)/2 | r > (B+1)*(1-alpha(a)/2)
            b_rejections = b_rejections + 1;
        end
    end
    fprintf("Percentage of rejections for H0: mX = mY at %.4f significance level: %.3f%%.\n", alpha(a), b_rejections/M*100);
    
    fprintf("RANDOM PERMUTATION\n")
    rp_rejections = 0;
    for j=1:M
        r = find(rp_dif_means(:,j) == muX(j)-muY(j));      %rank
        % if all the values are identical, select the middle rank
        if length(r) == B+1
            r = round((B+1)/2);
        elseif length(r) >= 2
            % If at least one bootstrap statistic is identical to the 
            % original, pick the rank of one of them at random
            r = r(unidrnd(length(r)));
        end
        if r < (B+1)*alpha(a)/2 | r > (B+1)*(1-alpha(a)/2)
            rp_rejections = rp_rejections + 1;
        end
    end
    fprintf("Percentage of rejections for H0: mX = mY at %.4f significance level: %.3f%%.\n", alpha(a), rp_rejections/M*100);
    fprintf("\n")
end

% Bootstrap hypothesis testing for mean difference
alpha = 0.05;
rej = 0;
B = 1000;       % Number of Bootstrap samples
h_bootstrap = zeros(M,1);
p_bootstrap = zeros(M,1);
for i=1:M
    b_dif_means = zeros(B+1,1);
    for j=1:B
        tmp_index = unidrnd(m+n, m+n, 1);
        tmp_data = [dataX(:,i); dataY(:,i)];
        bootstap_samplesX = tmp_data(tmp_index(1:n));
        bootstap_samplesY = tmp_data(tmp_index(n+1:m+n));
        b_dif_means(j) = mean(bootstap_samplesX)-mean(bootstap_samplesY);
    end
    b_dif_means(B+1) = muX(i)-muY(i);
    b_dif_means = sort(b_dif_means);
    r = find(b_dif_means == muX(i)-muY(i));      %rank
    % If all the values are identical, select the middle rank
    if length(r) == B+1
        r = round((B+1)/2);
    elseif length(r) >= 2
        % If at least one bootstrap statistic is identical to the 
        % original, pick the rank of one of them at random
        r = r(unidrnd(length(r)));
    end
    if r < (B+1)*alpha/2 | r > (B+1)*(1-alpha/2)
        h_bootstrap(i) = 1;
        rej = rej+1;
    else
        h_bootstrap(i) = 0;
    end
    if r > 0.5*(B+1)
        p_bootstrap(i) = 2*(1-r/(B+1));
    else
        p_bootstrap(i) = 2*r/(B+1);
    end
end