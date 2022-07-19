% Exercise 5.8

clc;
clear;
close all;

% Fit linear regression model using stepwise regression
physical_data = importdata('physical.txt');
data = physical_data.data;
var_names = char('Mass','Fore','Bicep','Chest','Neck','Shoulder', ...
    'Waist','Height','Calf','Thigh','Head');
Y = data(:,1);
X = data(:,2:end);
n = length(Y);
k = size(X,2);
alpha = 0.05;
zc = norminv(1-alpha/2);

% Full model using all the X variables
x = [ones(n,1) X];
[b,bint] = regress(Y,x);
% Predicted values
y = x * b;
% Error
e = Y - y;
se = sqrt((1/(n-(k+1)))*(sum(e.^2)));
estar = e/se;
% R-squared statistic
meanY = mean(Y);
R2 = 1-(sum(e.^2))/(sum((Y-meanY).^2));
adjR2 =1-((n-1)/(n-(k+1)))*(sum(e.^2))/(sum((Y-meanY).^2));

fprintf('Full model (10 variables)\n');
fprintf('Coefficients\n\tb0=%.5f, b0_int=[%.5f,%.5f]\n', ...
    b(1),bint(1,1),bint(1,2));
for i=2:length(b)
    fprintf('\t%s: b%d=%.5f, b%dint=[%.5f,%.5f]\n', ...
        var_names(i,:),i-1,b(i),i-1,bint(i,1),bint(i,2));
end
fprintf('Residual std: %.5f\n',se);
fprintf('R-squared: %.5f\nadj R-square: %.5f\n',R2,adjR2);

figure();
scatter(Y,estar,12,'blue','filled');
hold on;
yline(zc,'r-.');
yline(-zc,'r-.');
xlabel('Mass (kg)');
ylabel('Standardized adjustment error');
title('Diagnostic plot for full model');

% Stepwise regression model
[bs,sebs,~,finalmodel,stats]=stepwisefit(X,Y);
bs = [stats.intercept;bs];
model_vars = find(finalmodel==1);
% Predicted values
y = x * (bs.*transpose([1 finalmodel]));
% Error
e = Y - y;
k = sum(finalmodel);
se = sqrt((1/(n-(k+1)))*(sum(e.^2)));
estar = e/se;
% R-squared statistic
R2 = 1-(sum(e.^2))/(sum((Y-meanY).^2));
adjR2 =1-((n-1)/(n-(k+1)))*(sum(e.^2))/(sum((Y-meanY).^2));
tc = tinv(1-alpha/2,n-(k+1));
fprintf("Model using stepwise regression\n");
fprintf('Coefficients\n\tb0=%.5f\n',bs(1));
for i=1:k
    fprintf('\t%s: b%d=%.5f, b%dint=[%.5f,%.5f]\n', ...
        var_names(model_vars(i)+1,:),model_vars(i),bs(model_vars(i)+1), ...
        model_vars(i), bs(model_vars(i)+1)-tc*sebs(model_vars(i)+1), ...
        bs(model_vars(i)+1)+tc*sebs(model_vars(i)+1));
end
fprintf('Residual std: %.5f\n',se);
fprintf('R-squared: %.5f\nadj R-square: %.5f\n',R2,adjR2);

figure();
scatter(Y,estar,12,'blue','filled');
hold on;
yline(zc,'r-.');
yline(-zc,'r-.');
xlabel('Mass (kg)');
ylabel('Standardized adjustment error');
title('Diagnostic plot for stepwise regression model');