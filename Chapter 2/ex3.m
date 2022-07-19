% Exercise 2.3
clc;
clear;
% Prove that Var(X+Y) != Var(X) + Var(Y) when X and Y are not independent
muX=0;
muY=0;
sigma1=1;
sigma2=3;
corr=1.5;
n= 1000;

mu = [muX muY];
Sigma = [sigma1 corr; corr sigma2];
R = mvnrnd(mu,Sigma,n);

X = R(:,1);
Y = R(:,2); 

varX=var(X);
varY=var(Y);
varXsumY=var(X+Y);

fprintf("Var[X]+Var[Y]=%f and Var[X+Y]=%f.\n", varX+varY, varXsumY);