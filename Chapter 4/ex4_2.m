% Exercise 4.2

clc;
clear;
close all;

sigma = 5;      % uncertainty for length and width

% a
l = 500;
w = 300;
area_sigma = sqrt(w^2*sigma^2 + l^2*sigma^2);
fprintf("Uncertainty in area: %.4f\n", area_sigma);

l = [10:1000];
w = sqrt(area_sigma^2 - l.^2*sigma^2)/sigma;
figure();
plot(l,w);

% b
w = [10:1000];
n = length(w);
area_sigma = zeros(n);
for i = 1:n
    area_sigma(:,i) = sqrt(w.^2*sigma^2 + l(i)^2*sigma^2);
end
figure();
surf(l,w,area_sigma);
