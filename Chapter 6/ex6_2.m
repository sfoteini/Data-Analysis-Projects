% Exercise 6.2

clc;
clear;
close all;

x = importdata('yeast.dat.txt');
[n,p] = size(x);

figure();
plot(mean(x,2))
xlabel('Variable index i')
ylabel('Sample mean of x_i')
title('Sample mean of yeast data')
figure();
plot(std(x'));
xlabel('Variable index i')
ylabel('Sample SD of x_i')
title('Sample SD of yeast data')

% a : PCA and estimation of d < p
% Centering
y = x - repmat(sum(x)/n,n,1);
% Covariance matrix
covx = cov(y);
[eigvec,eigval] = eig(covx);
eigval = diag(eigval);          % The diagonal elements
% Sort in descending order
[eigval,ind] = sort(eigval,'descend');
eigvec = eigvec(:,ind);
% PC scores
z = y * eigvec;
% Scree plot
figure();
plot(1:p,eigval,'o-');
title('Scree Plot');
xlabel('Index');
ylabel('Eigenvalue');
title('Scree plot');
% Explained Variance percentage
td = 100*cumsum(eigval)/sum(eigval);
figure();
plot(1:p,td,'o-');
xlabel('Index');
ylabel('Variance percentage');
title('Explained variance percentage');
% Size of the variance
avgeig = mean(eigval);
ind = find(eigval > avgeig);
fprintf('Dimension d using size of the variance: %d\n',length(ind));
% b : Projection in 2-D and 3-D
z3 = y * eigvec(:,1:3);
figure();
scatter3(z3(:,1),z3(:,2),z3(:,3),10,'blue');
xlabel('PC1');
ylabel('PC2');
zlabel('PC3');
title('Principal component scores in 3-D');
z2 = y * eigvec(:,1:2);
figure();
scatter(z2(:,1),z2(:,2),10,'blue');
xlabel('PC1');
ylabel('PC2');
title('Principal component scores in 2-D');