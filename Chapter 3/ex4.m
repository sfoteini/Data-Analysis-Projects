% Exercise 3.4

clc;
clear;

data = [41 46 47 47 48 50 50 50 50 50 50 50 ...
        48 50 50 50 50 50 50 50 52 52 53 55 ...
        50 50 50 50 52 52 53 53 53 53 53 57 ...
        52 52 53 53 53 53 53 53 54 54 55 68];
oldVar = 25;
alpha = 0.05;
    
[hv,pv,civ] = vartest(data, oldVar, 'Alpha', alpha);
% a : Confidence Interval of variance
fprintf('%.2f %% Confidence Interval of variance: [%.4f, %.4f]\n', ...
    (1-alpha)*100, civ(1), civ(2));
fprintf('%.2f %% Confidence Interval of std: [%.4f, %.4f]\n', ...
    (1-alpha)*100, sqrt(civ(1)), sqrt(civ(2)));

% b : Null Hypothesis
if hv
    fprintf(['Rejection of the null Hypothesis H0: σ = %.2f at %.2f ' ...
        'significance level.\n'], sqrt(oldVar), alpha);
else
    fprintf(['Acceptance of the null Hypothesis H0: σ = %.2f at %.2f ' ...
        'significance level.\n'], sqrt(oldVar), alpha);
end
fprintf('p-value for H0: σ = %.2f: p = %.3f\n\n',sqrt(oldVar),pv);

% c : Confidence Interval of mean value
oldMean = 52;
[hm,pm,cim] = ttest(data, oldMean, 'Alpha', alpha);
fprintf('%.2f %% Confidence Interval of mean value: [%.4f, %.4f]\n', ...
    (1-alpha)*100, cim(1), cim(2));

% d : Null Hypothesis
if hm
    fprintf(['Rejection of the null Hypothesis H0: μ = %.2f at %.2f ' ...
        'significance level.\n'], oldMean, alpha);
else
    fprintf(['Acceptance of the null Hypothesis H0: μ = %.2f at %.2f ' ...
        'significance level.\n'], oldMean, alpha);
end
fprintf('p-value for H0: μ = %.2f: p = %.3f\n\n',oldMean,pm);

% e : Goodness-of-fit test and p-value
[h,p] = chi2gof(data, 'Alpha', alpha);
if h
    fprintf(['Rejection of the null Hypothesis H0: data comes from a ' ...
        'population with a normal distribution at %.2f significance ' ...
        'level.\n'], alpha);
else
    fprintf(['Acceptance of the null Hypothesis H0: data comes from a ' ...
        'population with a normal distribution at %.2f significance ' ...
        'level.\n'], alpha);
end
fprintf('The p-value of the goodness-of-fit test is %.3f.\n', p);