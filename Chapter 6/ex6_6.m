% Exercise 6.6

clc;
clear;
close all;

% Linear regression model using OLS, PCR, PLS, RR, LASSO, stepwise
% reression
physical_data = importdata('physical.txt');
data = physical_data.data;
var_names = char('Mass','Fore','Bicep','Chest','Neck','Shoulder', ...
    'Waist','Height','Calf','Thigh','Head');
y = data(:,1);
x = data(:,2:end);
n = length(y);
p = size(x,2);
alpha = 0.05;
zc = norminv(1-alpha/2);

% Centering the data
mux = mean(x);
xc = x - repmat(mux,n,1);
muy = mean(y);
yc = y - muy;

TSS = sum((y-muy).^2);
[u,sigma,v] = svd(xc,'econ');

% Ordinary Least Squares - OLS
bOLS = v * inv(sigma) * u' * yc;
bOLS = [muy - mux*bOLS; bOLS];
yOLS = [ones(n,1) x] * bOLS; 
resOLS = y-yOLS;
RSS_OLS = sum(resOLS.^2);
r2OLS = 1 - RSS_OLS/TSS;
% Scatter plot observed vs OLS calculated y
figure();
scatter(y,yOLS,8);
xlabel('y');
ylabel('OLS y');
title(sprintf('OLS Model R^2=%.4f',r2OLS));
% Standardized error
alpha = 0.05;
zc = norminv(1-alpha/2);
figure();
scatter(y,resOLS/std(resOLS),8);
hold on;
yline(zc,'r-.');
yline(-zc,'r-.');
xlabel('y');
ylabel('Standardized error');
title('OLS standardized errors');

% Principal Component Regression - PCR
d = 2;                  % dimension reduction
lambda = zeros(p,1);
lambda(1:d) = 1;
bPCR = v * diag(lambda) * inv(sigma) * u'* yc;
bPCR = [muy - mux*bPCR; bPCR];
yPCR = [ones(n,1) x] * bPCR;
resPCR = y - yPCR;
RSS_PCR = sum(resPCR.^2);
r2PCR = 1 - RSS_PCR/TSS;
figure();
scatter(y,yPCR,8);
xlabel('y');
ylabel('PCR y');
title(sprintf('PCR Model R^2=%.4f',r2PCR));
figure();
scatter(y,resPCR/std(resPCR),8);
hold on
yline(zc,'r-.');
yline(-zc,'r-.');
xlabel('y');
ylabel('Standardized error');
title('PCR standardized errors');

% Principal Least Squares - PLS
[Xl,Yl,Xscores,Yscores,bPLS] = plsregress(x,y,d);
yPLS = [ones(n,1) x]*bPLS;
resPLS = y - yPLS;
RSS_PLS = sum(resPLS.^2);
r2PLS = 1 - RSS_PLS/TSS;
figure();
scatter(y,yPLS,8);
xlabel('y');
ylabel('PLS y');
title(sprintf('PLS Model R^2=%.4f',r2PLS));
figure();
scatter(y,resPLS/std(resPLS),8);
hold on;
yline(zc,'r-.');
yline(-zc,'r-.');
xlabel('y');
ylabel('Standardized error');
title('PLS standardized errors');

% Ridge regression
mu = RSS_OLS/(n-p);
sig = diag(sigma);
lambda = sig.^2 ./ (sig.^2 + mu);
bRR = v * diag(lambda) * inv(sigma) * u'* yc;
bRR = [muy - mux*bRR; bRR];
yRR = [ones(n,1) x] * bRR; 
resRR = y - yRR;
RSS_RR = sum(resRR.^2);
r2RR = 1 - RSS_RR/TSS;
figure();
scatter(y,yRR,8);
xlabel('y');
ylabel('RR y');
title(sprintf('RR Model R^2=%.4f',r2RR));
figure();
scatter(y,resRR/std(resRR),8);
hold on
yline(zc,'r-.');
yline(-zc,'r-.');
xlabel('y');
ylabel('Standardized error');
title('RR standardized errors');

% LASSO
[bL,fitinfo] = lasso(xc,yc);
lassoPlot(bL,fitinfo,'PlotType','Lambda','XScale','log');
lambda = 0.5;
[lmin, ilmin] = min(abs(fitinfo.Lambda - lambda));
bLASSO = bL(:,ilmin);
bLASSO = [muy - mux*bLASSO; bLASSO];
yLASSO = [ones(n,1) x] * bLASSO;
resLASSO = y - yLASSO;
RSS_LASSO = sum(resLASSO.^2);
r2LASSO = 1 - RSS_LASSO/TSS;
figure();
scatter(y,yLASSO,8);
hold on;
xlabel('y')
ylabel('LASSO y')
title(sprintf('LASSO Model R^2=%.4f',r2LASSO));
figure();
scatter(y,resLASSO/std(resLASSO),8);
hold on;
yline(zc,'r-.');
yline(-zc,'r-.');
xlabel('y');
ylabel('Standardized error');
title('LASSO standardized errors');

% Stepwise regression model
[bs,sebs,~,finalmodel,stats]=stepwisefit(x,y,'Display','off');
bs = [stats.intercept;bs];
model_vars = find(finalmodel==1);
bSTEP = bs.*transpose([1 finalmodel]);
ySTEP = [ones(n,1) x] * bSTEP;
resSTEP = y-ySTEP;
RSS_STEP = sum(resSTEP.^2);
r2STEP = 1 - RSS_STEP/TSS;
figure();
scatter(y,ySTEP,8);
xlabel('y');
ylabel('Stepwise y');
title(sprintf('Stepwise Model R^2=%.4f',r2STEP));
figure();
scatter(y,resSTEP/std(resSTEP),8);
hold on;
yline(zc,'r-.');
yline(-zc,'r-.');
xlabel('y');
ylabel('Standardized error');
title('Stepwise standardized errors');

% Coefficients
fprintf('\t\t OLS \t\t PCR \t\t PLS \t\t RR \t\t LASSO \t\t STEP\n');
fprintf('Constant \t %5.3f \t %5.3f \t %5.3f \t %5.3f \t %5.3f \t %5.3f \n',...
    bOLS(1),bPCR(1),bPLS(1),bRR(1),bLASSO(1),bSTEP(1));
for i=2:p+1
    fprintf('%s \t %5.3f \t\t %5.3f \t\t %5.3f \t\t %5.3f \t\t %5.3f \t\t %5.3f\n',...
        var_names(i,:),bOLS(i),bPCR(i),bPLS(i),bRR(i),bLASSO(i),bSTEP(i));
end