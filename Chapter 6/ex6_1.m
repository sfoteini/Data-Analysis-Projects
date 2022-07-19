% Exercise 6.1

clc;
clear;
close all;

% Generate n random numbers from the same multivariate normal distribution.
n = 1000;
mu = [0 0];
sigma = [1 0; 0 4];
x = mvnrnd(mu', sigma, n);
w = [0.2 0.8; 0.4 0.5; 0.7 0.3];
x3 = x*w';       % 3-D data points
p = size(w,1);

% Plot the data points in 2-D and 3-D
figure();
scatter(x(:,1),x(:,2),10,'blue');
xlabel('x1');
ylabel('x2');
title('2D Gaussian generated points');
figure();
scatter3(x3(:,1),x3(:,2),x3(:,3),10,'blue');
xlabel('y1');
ylabel('y2');
ylabel('y3');
title('3D observed points');

% a : PCA (Eigenvalues and eigenvectors of the covariance matrix and
% principal components scores)
% Centering
y = x3 - repmat(sum(x3)/n,n,1);
% Covariance matrix
covx = cov(y);
[eigvec,eigval] = eig(covx);
eigval = diag(eigval);          % The diagonal elements
% Sort in descending order
[eigval,ind] = sort(eigval,'descend');
eigvec = eigvec(:,ind);
% PC scores
z3 = (eigvec * y')';

% b : Scree plot
figure();
plot(1:p,eigval,'o-');
title('Scree Plot');
xlabel('Index');
ylabel('Eigenvalue');
title('Scree plot');
% Plot PC scores
figure();
scatter3(z3(:,1),z3(:,2),z3(:,3),10,'blue');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('Principal component scores');

% c : PC scores in 2D
z2 = (eigvec(1:2,:) * y')';
figure();
scatter(z2(:,1),z2(:,2),10,'blue');
xlabel('PC1');
ylabel('PC2');
title('Principal component scores');