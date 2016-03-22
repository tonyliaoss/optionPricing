function [ C ] = lseRegression( X, Y )
% lseRegression: performs the least-squares-error regression
% with the set of basis functions [1, x, x*2]
%   Return values: a*x^2 + b*x + c

if isempty(X)
  C = [0 0 0];
else
  C = polyfit(X, Y, 2);
end

end
