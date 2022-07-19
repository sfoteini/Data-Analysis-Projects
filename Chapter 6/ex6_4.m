% Exercise 6.4

clc;
clear;
close all;

% Mixing of two sound signals and decomposition with ICA
n = 1e4;                             % number of samples
s1 = load('chirp.mat');
s1 = s1.y;
s1 = s1(1:n);
s2 = load('gong.mat');
s2 = s2.y;
s2 = s2(1:n);
s = [s1 s2];
s(:,1) = s(:,1) - mean(s(:,1));
s(:,2) = s(:,2) - mean(s(:,2));
% Plot the 2 signals
figure();
plot(s(:,1));
xlabel('t');
ylabel('s1(t)');
title('Source signal 1: chirp');
figure();
plot(s(:,2));
xlabel('t');
ylabel('s2(t)');
title('Source signal 2: gong');
% Scatter
figure();
scatter(s(:,1),s(:,2),8);
xlabel('s1');
ylabel('s2');
title("Scatter plot of source signals");

% Mix the signals and apply ICA with and without prewhitening
prewhitening = 1;
% Mixing matrix A
% a : matrix size 2x2
A = [-0.1 0.3; -2.5 -0.2];
% b : matrix size 2x3
%A = [-0.1 0.3 -0.2; -2.5 -0.2 -0.3];
p = size(A,2);
x = s * A;
% Plot the mixed signals
figure();
plot(x(:,1));
xlabel('t');
ylabel('x1(t)');
title('Mixed signal 1');
figure();
plot(x(:,2));
xlabel('t');
ylabel('x2(t)');
title('Mixed signal 2');
if p == 3
    figure();
    plot(x(:,3));
    xlabel('t');
    ylabel('x3(t)');
    title('Mixed signal 3');
end
% Centering
y = x - repmat(sum(x)/n,n,1);
% Scatter
if p == 2
    figure();
    scatter(y(:,1),y(:,2),8);
    xlabel('y1');
    ylabel('y2');
    title("Scatter plot of observed signals");
elseif p == 3
    figure();
    scatter3(y(:,1),y(:,2),y(:,3),8);
    xlabel('y1');
    ylabel('y2');
    zlabel('y3');
    title("Scatter plot of observed signals");
end

% Prewhitening
if prewhitening
    y = prewhiten(y);
end

% ICA
Mdl = rica(y,p,'Lambda',0.5);
w = Mdl.TransformWeights;
z = transform(Mdl,y);
% Plot the reconstructed signals
figure();
plot(z(:,1));
xlabel('t');
ylabel('s1(t)');
title('ICA reconstructed s1');
figure();
plot(z(:,2));
xlabel('t');
ylabel('s2(t)');
title('ICA reconstructed s2');
% Scatter
if p == 2
    figure();
    scatter(z(:,1),z(:,2),8);
    xlabel('s1');
    ylabel('s2');
    title("Scatter plot of ICA reconstructed signals");
elseif p == 3
    figure();
    plot(z(:,3));
    xlabel('t');
    ylabel('s3(t)');
    title('ICA reconstructed s3');
    figure();
    scatter3(z(:,1),z(:,2),z(:,3),8);
    xlabel('s1');
    ylabel('s2');
    zlabel('s3');
    title("Scatter plot of ICA reconstructed signals");
end


% source: https://www.mathworks.com/help/stats/extract-mixed-signals.html
function y = prewhiten(x)
    % Size of x
    [n,p] = size(x);
    assert(n >= p);
    % SVD of covariance of x
    [U,Sig] = svd(cov(x));
    Sig = diag(Sig);
    Sig = Sig(:)';
    % Figure out which values of Sig are non-zero
    tol = eps(class(x));
    idx = (Sig > max(Sig)*tol);
    assert(~all(idx == 0));
    % Get the non-zero elements of Sig and corresponding columns of U
    Sig = Sig(idx);
    U = U(:,idx);
    % Compute prewhitened data
    mu = mean(x,1);
    y = bsxfun(@minus,x,mu);
    y = bsxfun(@times,y*U,1./sqrt(Sig));
end