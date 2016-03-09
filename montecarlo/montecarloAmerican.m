function [ P ] = montecarloAmerican( S, tau, E, r, sigma )
% Gives the price of an American put option without dividends, as computed
% through Monte-Carlo simulation by least squares (Longstaff and Schwartz).
%   Detailed explanation goes here
% INPUT PARAMETERS
%   S - the current market value of the asset.
%   tau - time to expiry of the option, expressed in years.
%   E - the strike price, or exercise price, of the option.
%   r - the risk-free investment return per annum.
%   sigma - the market volatility.


d1 = (log(S/E) + (r + 0.5 * sigma ^ 2) * tau) / (sigma * sqrt(tau));

d2 = (log(S/E) + (r - 0.5 * sigma ^ 2) * tau) / (sigma * sqrt(tau));

P = E * exp(-r * tau) * normcdf(-d2) - S * normcdf(-d1);

%% note also need to choose number and type of basis functions.
% e.g. number = 2 and type = polynomial -> least squares.

% default values for optional arguments...
defaultN = 2000; % number of sample paths.
defaultNumPartitionsT = 2000; % number of observation dates.
defaultZeta = 1e6; % zeta is the penalty parameter.

% ballpark: we need to draw (N * numPartitionsT) normally distributed random
% numbers.

end

