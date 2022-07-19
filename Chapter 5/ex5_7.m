% Exercise 5.7

clc;
clear;
close all;

% Data
R = [0.76 0.86 0.97 1.11 1.45 1.67 1.92 2.23 2.59 3.02 3.54 4.16 4.91 ...
    5.83 6.94 8.31 10.00 12.09 14.68 17.96 22.05 27.28 33.89 42.45 ...
    53.39 67.74 86.39 111.30 144.00 188.40 247.50 329.20]';
T = [110 105 100 95 85 80 75 70 65 60 55 50 45 40 35 30 25 20 15 10 5 ...
    0 -5 -10 -15 -20 -25 -30 -35 -40 -45 -50]' + 273.15;
lnR = log(R);
Tinverse = 1./T;
n = length(R);
alpha = 0.05;

% a : Find the polynomial regression model
kmax = 4;           % max degree of the polynomial regression model
% Correlation coefficient
r = corrcoef(lnR,Tinverse);
r = r(1,2);
fprintf('Correlation coefficient: %.3f\n',r);

% Diagnostic plot
npoints = 100;
R0 = linspace(min(lnR),max(lnR),npoints)';
zc = norminv(1-alpha/2);
for k=1:kmax
    % Choose a polynomial model
    switch k
        case 1
            x = [ones(n,1) lnR];
            x0 = [ones(npoints,1) R0];
        case 2
            x = [ones(n,1) lnR lnR.^2];
            x0 = [ones(npoints,1) R0 R0.^2];
        case 3
            x = [ones(n,1) lnR lnR.^2 lnR.^3];
            x0 = [ones(npoints,1) R0 R0.^2 R0.^3];
        case 4
            x = [ones(n,1) lnR lnR.^2 lnR.^3 lnR.^4];
            x0 = [ones(npoints,1) R0 R0.^2 R0.^3 R0.^4];
    end
    % Fit the regression model
    [b, bint,r,rint,stats] = regress(Tinverse,x);
    % Predicted values of 1/T
    y = x * b;
    y0 = x0 * b;
    % Error
    e = Tinverse - y;
    se = sqrt((1/(n-(k+1)))*(sum(e.^2)));
    estar = e./se;
    % R-squared statistic
    meanY = mean(Tinverse);
    R2 = 1-(sum(e.^2))/(sum((Tinverse-meanY).^2));
    adjR2 =1-((n-1)/(n-(k+1)))*(sum(e.^2))/(sum((Tinverse-meanY).^2));
    % Linear model and diagnostic plot
    figure();
    subplot(2,1,1);
    scatter(lnR,Tinverse,12,'blue','filled');
    title(sprintf('Ln of R vs 1/T: Degree = %d', k));
    hold on;
    plot(R0,y0,'LineWidth',2);
    xlabel('Ln of R');
    ylabel('1/T');
    txt = ['R2=' num2str(R2) newline 'djR2=' num2str(adjR2)];
    text(R0(round(npoints/2)),y0(round(npoints/5)),txt)
    hold off;
    subplot(2,1,2);
    scatter(Tinverse,estar,12,'blue','filled');
    hold on;
    yline(zc,'r-.');
    yline(-zc,'r-.');
    xlabel('1/T');
    ylabel('Standardized adjustment error');
    title('Diagnostic plot');
    % Print the results
    fprintf('Polynomial degree=%d\n',k);
    for i=1:k+1
        fprintf('\tb_%d = %.10f\n',i-1,b(i));
    end
end

% b : Steinhart-Hart's model
x = [ones(n,1) lnR lnR.^3];
x0 = [ones(npoints,1) R0 R0.^3];
% Fit the regression model
[b, bint,r,rint,stats] = regress(Tinverse,x);
% Predicted values of 1/T
y = x * b;
y0 = x0 * b;
% Error
e = Tinverse - y;
se = sqrt((1/(n-2))*(sum(e.^2)));
estar = e./se;
% R-squared statistic
R2 = 1-(sum(e.^2))/(sum((Tinverse-meanY).^2));
adjR2 =1-((n-1)/(n-(k+1)))*(sum(e.^2))/(sum((Tinverse-meanY).^2));

figure();
subplot(2,1,1);
scatter(lnR,Tinverse,12,'blue','filled');
title('Steinhart-Hart model');
hold on;
plot(R0,y0,'LineWidth',2);
xlabel('Ln of R');
ylabel('1/T');
txt = ['R2=' num2str(R2) newline 'djR2=' num2str(adjR2)];
text(R0(round(npoints/2)),y0(round(npoints/5)),txt)
hold off;
subplot(2,1,2);
scatter(Tinverse,estar,12,'blue','filled');
hold on;
yline(zc,'r-.');
yline(-zc,'r-.');
xlabel('1/T');
ylabel('Standardized adjustment error');
title('Diagnostic plot');
% Print the results
fprintf('Steinhart-Hart model\n');
for i=1:3
    fprintf('\tb_%d = %.10f\n',i-1,b(i));
end