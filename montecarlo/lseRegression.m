function [ C ] = lseRegression( X, Y )
% lseRegression: performs the least-squares-error regression
% with the set of basis functions [1, x, x*2]
%   Return values: a*x^2 + b*x + c

% transform the data before polyfitting.
[p, S, mu] = polyfit(X, Y, 2);
% re-transform the data.
aa = p(1);
bb = p(2);
cc = p(3);
mu1 = mu(1);
mu2 = mu(2);

c = aa * mu1^2/mu2^2 - (bb * (mu1 / mu2)) + cc;
b = (bb * mu2 - 2 * aa * mu1) / (mu2^2);
a = aa / mu2 ^2;

C = [a b c];

end
