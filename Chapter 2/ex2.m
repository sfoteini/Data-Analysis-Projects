% Exercise 2.2
clc;
clear;
close all;
% Generate random numbers from exponential distribution 
lamda = 1;
uniform_rnd = rand(1e3, 1);
exp_rnd = -(1/lamda).*log(1-uniform_rnd);

x = 0.01:0.1:10;
exp_pdf = lamda.*exp(-lamda.*x);

histogram(exp_rnd, 'Normalization','pdf');
hold on
plot(x, exp_pdf, 'LineWidth',1.5);
legend('simulated','analytic');
title('Exponential distribution from random numbers');