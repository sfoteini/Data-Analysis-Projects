% Exercise 3.3

clc;
clear;
% a
expMeanConfInterval(1/15, 1000, 5);
% b
expMeanConfInterval(1/15, 1000, 100);

function expMeanConfInterval(lambda, M, n)
   meanX = 1/lambda;
   rejection = 0;
   X = exprnd(meanX, n, M);
   % Confidence Interval (a=0.05)
   ci = zeros(2,M);
   for i=1:M
    [h,p,ci(:,i)] = ttest(X(:, i),meanX);
    %fprintf("Confidence Interval: [%f, %f]\n", ci(1), ci(2));
    if h
        rejection = rejection + 1;
    end
   end
   fprintf("For M=%d samples and sample size n=%d the mean = %.4f is " + ...
       "inside confidence interval with percentage %.4f %%.\n", ...
       M, n, meanX, (1-rejection/M)*100);
   figure();
   histogram(ci(1,:));
   title('Histogram of lower limit of ci');
   figure();
   histogram(ci(2,:));
   title('Histogram of upper limit of ci');
end