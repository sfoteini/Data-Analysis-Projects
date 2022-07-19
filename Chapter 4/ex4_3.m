% Exercise 4.3

clc;
clear;
close all;

M = 1000;
muv = 77.78;
sigmav = 0.71;
mui = 1.21;
sigmai = 0.071;
muf = 0.283;
sigmaf = 0.017;
rvf = 0.5;

% a
sigmap = sqrt((mui*cos(muf))^2*sigmav^2 + (muv*cos(muf))^2*sigmai^2 + ...
    (muv*mui*(-sin(muf)))^2*sigmaf^2);
fprintf("Uncertainty of P: %.4f\n", sigmap);

% b
V = normrnd(muv, sigmav, M, 1);
I = normrnd(mui, sigmai, M, 1);
f = normrnd(muf, sigmaf, M, 1);
P = V.*I.*cos(f);
stdp = std(P);
fprintf("Uncertainty of P based on the experiments: %.4f\n", stdp);

% c
covvf = rvf*sigmav*sigmaf;
sigmap2 = sqrt((mui*cos(muf))^2*sigmav^2 + (muv*cos(muf))^2*sigmai^2 + ...
    (muv*mui*(-sin(muf)))^2*sigmaf^2 + 2*mui*cos(muf)*muv*mui*(-sin(muf))*covvf);
sigma = [sigmav^2 0 covvf; 0 sigmai^2 0; covvf 0 sigmaf^2];
mu = [muv mui muf];
data = mvnrnd(mu, sigma, M);
V2 = data(:, 1);
I2 = data(:, 2);
f2 = data(:, 3);
P2 = V2.*I2.*cos(f2);
stdp2 = std(P2);
fprintf("Uncertainty of P: %.4f\n", sigmap2);
fprintf("Uncertainty of P based on the experiments: %.4f\n", stdp2);