function [ P ] = BSEqnEuropean( S, tau, E, r, sigma )
% BsEqnEuropean: Gives the price of a European put option without
% dividends, as computed by the closed-form Black-Scholes equations.
%
% INPUT PARAMETERS
%   S - the current market value of the asset.
%   tau - time to expiry of the option, expressed in years.
%   E - the strike price, or exercise price, of the option.
%   r - the risk-free investment return per annum.
%   sigma - the market volatility.

d1 = (log(S/E) + (r + 0.5 * sigma ^ 2) * tau) / (sigma * sqrt(tau));

d2 = (log(S/E) + (r - 0.5 * sigma ^ 2) * tau) / (sigma * sqrt(tau));

P = E * exp(-r * tau) * normcdf(-d2) - S * normcdf(-d1);

end

