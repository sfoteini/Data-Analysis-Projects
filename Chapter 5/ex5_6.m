% Exercise 5.6

clc;
clear;
close all;

% Data
X = [2 3 8 16 32 48 64 80]';                        % distance
Y = [98.2 91.7 81.3 64.0 36.4 32.6 17.1 11.3]';     % percentage
n = length(X);
alpha = 0.05;

% a : scatter plot, linear regression model and diagnostic plot
r = corrcoef(X,Y);
r = r(1,2);
fprintf('Correlation coefficient: %.3f\n',r);

b = regress(Y,[ones(n,1) X]);
fprintf('Linear regression model: y = ax + b -> Y = %.4f X + %.4f\n', ...
    b(2),b(1));
x0 = linspace(min(X),max(X),100)';
y0 = b(2)*x0 + b(1);
figure();
scatter(X,Y,12,'blue','filled');
title(sprintf('Km driven vs Percentage usable: r=%.3f', r));
hold on;
plot(x0,y0,'LineWidth',2);
xlabel('Km driven');
ylabel('Percentage usable');

e = Y - [ones(n,1) X] * b;
se = sqrt((1/(n-2))*(sum(e.^2)));
estar = e./se;
zc = norminv(1-alpha/2);
figure();
scatter(Y,estar,12,'blue','filled');
hold on;
yline(zc,'r-.');
yline(-zc,'r-.');
xlabel('Percentage usable');
ylabel('Standardized adjustment error');
title('Diagnostic plot');

% y = aexp(-bx)
lnY = log(Y);
varlnY = var(lnY);
mlnY = mean(lnY);
varX = var(X);
[bln,bintln] = regress(lnY,[ones(n,1) X]);
seln = sqrt((n-1)/(n-2)*(varlnY-bln(2)^2*varX));
eln = lnY - [ones(n,1) X] * bln;
R2ln = 1-(sum(eln.^2))/(sum((lnY-mlnY).^2));
adjR2ln = 1-((n-1)/(n-2))*(sum(eln.^2))/(sum((lnY-mlnY).^2));
fprintf(['Ln-transformed linear regression model ' ...
    'y = ln(a) + bx -> y = %.3f + (%.3f)x\n'],bln(1),bln(2));
fprintf('se = %.4f, R^2 = %.4f, adjR^2 = %.4f\n',seln,R2ln,adjR2ln);
lny0 = bln(2)*x0 + bln(1);
figure();
scatter(X,lnY,12,'blue','filled');
title(sprintf('Km driven vs ln of Percentage usable y=%.3f+(%.3f)x', ...
    bln(1),bln(2)));
hold on;
plot(x0,lny0,'LineWidth',2);
xlabel('Km driven');
ylabel('ln of Percentage usable');

% Diagnostic plot
estarln = eln./seln;
figure();
scatter(lnY,estarln,12,'blue','filled');
hold on;
yline(zc,'r-.');
yline(-zc,'r-.');
xlabel('ln of Percentage usable');
ylabel('Standardized adjustment error');
title('Diagnostic plot of ln linear regession');

% b : Percentage usable Y for X = 25
a = exp(bln(1));
b = bln(2);
fprintf('Exponential model: y = %.3f exp(%.3fx)\n',a,b);
y0 = a*exp(b*x0);
figure();
scatter(X,Y,12,'blue','filled');
title('Km driven vs Percentage usable');
hold on;
plot(x0,y0,'LineWidth',2);
xlabel('Km driven');
ylabel('Percentage usable');
text(mean(X),mean(Y),sprintf('y = %.3f exp(%.3fx)',a,b));
x = 25;
y = a*exp(b*x);
fprintf('Percentage usable for x = %.1fKm: y = %.4f',x*1000,y);