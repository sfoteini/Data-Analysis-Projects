% Exercise 5.4

clc;
clear;
close all;

air_data = importdata("lightair.dat.txt");
X = air_data(:,1);      % density
Y = air_data(:,2);      % light speed
n = length(X);
alpha = 0.05;
% Data exploration
mX = mean(X);
mY = mean(Y);
varX = var(X);
varY = var(Y);
covXY = cov(X,Y);
covXY = covXY(1,2);
fprintf("Air density:\n\tMean: %.3f\n\tVariance: %.3f\n",mX,varX);
fprintf("Light speed:\n\tMean: %.3f\n\tVariance: %.3f\n",mY,varY);
fprintf("Covariance: %.3f\n",covXY);
subplot(2,1,1);
histogram(X);
title("Air Density");
subplot(2,1,2);
histogram(Y);
title("Light speed");

% a : Correlation coeficient and scatter plot
r = corrcoef(X,Y);
r = r(1,2);
fprintf("Correlation coefficient: %.3f\n",r);
figure();
scatter(X,Y,12,"blue","filled");
title(sprintf("Air density and light speed r=%.3f", r));
xlabel("Air Density");
ylabel("Light speed");

% b : Linear regression model using least squares method
[b,bint] = regress(Y,[ones(n,1) X]);
% Least squares error
se = sqrt((n-1)/(n-2)*(varY-b(2)^2*varX));
fprintf("Linear regression model: y = ax + b -> " + ...
    "Y = %.4f X + %.4f\n",b(2),b(1));
fprintf("Confidence interval for intercept = %.4f: " + ...
    "[%.4f, %.4f]\n",b(1),bint(1,1),bint(1,2));
fprintf("Confidence interval for regression coefficient = %.4f: " + ...
    "[%.4f, %.4f]\n",b(2),bint(2,1),bint(2,2));
fprintf("Least squares error: se = %.4f\n",se);

% c : Draw the least squares line on a scatter plot
x0 = 1.29;      % air density
y0 = b(2)*x0 + b(1);
tc = tinv(1-alpha/2,n-2);
% Confidence interval for the mean value of Y for X=x0
Sxx = varX*(n-1);
mystd = se * sqrt(1/n + (x0-mX)^2/Sxx);
myci = [y0-tc*mystd y0+tc*mystd];
fprintf("Confidence interval for the mean value of Y for X=%.4f: " + ...
    "[%.4f, %.4f]\n",x0,myci(1),myci(2));
% Forecast interval for Y=y0, X=x0
y0std = se*sqrt(1+1/n+(x0-mX)^2/Sxx);
y0ci = [y0-tc*y0std y0+tc*y0std];
fprintf("Forecast interval for Y=%.4f and X=%.4f: " + ...
    "[%.4f, %.4f]\n",y0,x0,y0ci(1),y0ci(2));
% Calculate ci for every (x,y) and display them on the scatter plot
x0 = linspace(min(X),max(X),100)';
y0 = b(2)*x0 + b(1);
mystd = se * sqrt(1/n + (x0-mX).^2/Sxx);
myci = [y0-tc*mystd  y0+tc*mystd];
y0std = se*sqrt(1+1/n + (x0-mX).^2/Sxx);
y0ci = [y0-tc*y0std  y0+tc*y0std];

figure();
scatter(X,Y,10,"blue","filled");
hold on;
plot(x0,y0,'LineWidth',2);
xlabel("Air Density");
ylabel("Light speed");
if b(2)<0
    txt = sprintf("y = %.3f - %.3fx",b(1),abs(b(2))); 
else
    txt = sprintf("y = %.3f + %.3fx",b(1),abs(b(2))); 
end
title(strcat("Air density and light speed: ", txt));
plot(x0,y0ci(:,1),'m--');
plot(x0,y0ci(:,2),'m--');
plot(x0,myci(:,1),'c--');
plot(x0,myci(:,2),'c--');
legend("","Least squares line","Forecast interval","", ...
    "Mean Y confintence  interval");

% d : Real linear model
c = [299792.458-299000 -299792.458*0.00029/1.29];
fprintf("Real linear model: y = %.4fx + %.4f\n",c(2),c(1));
if c(1) < bint(1,1) | c(1) > bint(1,2)
    fprintf("Real intercept %.4f is not on the confidence interval " + ...
        "[%.4f, %.4f].\n",c(1),bint(1,1),bint(1,2));
else
    fprintf("Real intercept %.4f is on the confidence interval " + ...
        "[%.4f, %.4f].\n",c(1),bint(1,1),bint(1,2));
end
if c(2) < bint(2,1) | c(2) > bint(2,2)
    fprintf("Real slope %.4f is not on the confidence interval " + ...
        "[%.4f, %.4f].\n",c(2),bint(2,1),bint(2,2));
else
    fprintf("Real slope %.4f is on the confidence interval " + ...
        "[%.4f, %.4f].\n",c(2),bint(2,1),bint(2,2));
end
% light speed using the real model
v = c(1) + c(2).*x0;
accept = zeros(length(x0),1);
for i = 1:length(x0)
    if v(i) >= myci(i,1) & v(i) <= myci(i,2)
        accept(i) = 1;
    end
end
fprintf("The real light speed is in the confidence interval: " + ...
    "%.4f %%\n",sum(accept)*100/length(x0));