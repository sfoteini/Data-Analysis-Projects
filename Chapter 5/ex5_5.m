% Exercise 5.5

clc;
clear;
close all;

air_data = importdata("lightair.dat.txt");
X = air_data(:,1);      % density
Y = air_data(:,2);      % light speed
n = length(X);
alpha = 0.05;
M = 1000;

% Linear regression model using least squares method
[b,bint] = regress(Y,[ones(n,1) X]);

% Bootstrap confidence interval for regression coefficients
b_boot = zeros(M,2);
for i = 1:M
    index = unidrnd(n,n,1);
    Xb = X(index);
    Yb = Y(index);
    b_boot(i,:) = regress(Yb,[ones(n,1) Xb])';
end
% Lower and upper limit for ci
lowlim = round(M*alpha/2);
upplim = round(M*(1-alpha/2));
% Sort the b_boot array
b_boot_sorted = sort(b_boot);
bci = [b_boot_sorted(lowlim,1) b_boot_sorted(upplim,1);
       b_boot_sorted(lowlim,2) b_boot_sorted(upplim,2)];

fprintf("Linear model: intercept b0 = %.3f, slope b1 = %.3f\n",b(1),b(2));
fprintf("Parametric confintece interval for intercept: " + ...
    "[%.3f, %.3f]\n",bint(1,1),bint(1,2));
fprintf("Bootstrap confintece interval for intercept: " + ...
    "[%.3f, %.3f]\n",bci(1,1),bci(1,2));
fprintf("Parametric confintece interval for slope: " + ...
    "[%.3f, %.3f]\n",bint(2,1),bint(2,2));
fprintf("Bootstrap confintece interval for slope: " + ...
    "[%.3f, %.3f]\n",bci(2,1),bci(2,2));