% Exercise 5.9

clc;
clear;
close all;

% Fit linear regression model using stepwise regression
hospital_data = importdata('hospital.txt');
data = hospital_data.data;
var_names = char('ManHours','Cases','Eligible','OpRooms');
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

fprintf('Full model (3 variables)\n');
fprintf('Coefficients\n\tb0=%.5f, b0_int=[%.5f,%.5f]\n', ...
    b(1),bint(1,1),bint(1,2));
for i=2:length(b)
    fprintf('\t%s: b%d=%.5f, b%dint=[%.5f,%.5f]\n', ...
        var_names(i,:),i-1,b(i),i-1,bint(i,1),bint(i,2));
end
fprintf('Residual std: %.5f\n',se);
fprintf('R-squared: %.5f\nadj R-square: %.5f\n\n',R2,adjR2);

figure();
scatter(Y,estar,12,'blue','filled');
hold on;
yline(zc,'r-.');
yline(-zc,'r-.');
xlabel('Man Hours');
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
k1 = sum(finalmodel);
se = sqrt((1/(n-(k1+1)))*(sum(e.^2)));
estar = e/se;
% R-squared statistic
R2 = 1-(sum(e.^2))/(sum((Y-meanY).^2));
adjR2 =1-((n-1)/(n-(k1+1)))*(sum(e.^2))/(sum((Y-meanY).^2));
tc = tinv(1-alpha/2,n-(k1+1));
fprintf('Model using stepwise regression\n');
fprintf('Coefficients\n\tb0=%.5f\n',bs(1));
for i=1:k1
    fprintf('\t%s: b%d=%.5f, b%dint=[%.5f,%.5f]\n', ...
        var_names(model_vars(i)+1,:),model_vars(i),bs(i),model_vars(i), ...
        bs(model_vars(i))-tc*sebs(model_vars(i)), ...
        bs(model_vars(i))+tc*sebs(model_vars(i)));
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

% Check for multicollinearity
for i=1:k
    vars = setdiff(1:k,i);
    x = [ones(n,1) X(:,vars)];
    [b,bint] = regress(X(:,i),x);
    % Predicted values
    y = x * b;
    % Error
    e = X(:,i) - y;
    se = sqrt((1/(n-(k+1)))*(sum(e.^2)));
    estar = e/se;
    % R-squared statistic
    meanXi = mean(X(:,i));
    R2 = 1-(sum(e.^2))/(sum((X(:,i)-meanXi).^2));
    adjR2 =1-((n-1)/(n-(k+1)))*(sum(e.^2))/(sum((X(:,i)-meanXi).^2));

    fprintf('\nRegression model of %s\n',var_names(i+1,:));
    fprintf('\tResidual std: %.5f\n',se);
    fprintf('\tR-squared: %.5f\n\tadj R-square: %.5f\n',R2,adjR2);
end
figure();
plotmatrix(X);