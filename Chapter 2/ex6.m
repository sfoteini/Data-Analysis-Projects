% Exercise 2.6

clc;
clear;
close all;
% Simulation of the central limit theorem
nOfVars = 100;
nOfSamples = 1e4;
X_uniform = rand(nOfVars, nOfSamples);

Y = mean(X_uniform);

h = histfit(Y);
legend("Y distribution", "Normal distribution")