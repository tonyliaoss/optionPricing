function [ S ] = simulateBrownian( S0, mu, sigma, dt, nsteps, nsims )
% Function to generate sample paths for assets assuming geometric
% Brownian motion.
%
% S = simulateBrownian(S0,mu,sig,dt,steps,nsims)
%
% Inputs: S0 - stock price
%       : mu - expected return
%       : sig - volatility
%       : dt - size of time steps
%       : steps - number of time steps to calculate
%       : nsims - number of simulation paths to generate
%
% Output: S - a matrix where each column represents a simulated
%             asset price path.
%
% Notes: This code focuses on details of the implementation of the
%        Monte-Carlo algorithm.
%        It does not contain any programatic essentials such as error
%        checking.
%        It does not allow for optional/default input arguments.
%        It is not optimized for memory efficiency or speed.

% Author: Phil Goddard (phil@goddardconsulting.ca)
% Date: Q2, 2006

% Modified by Dazhi Liao (t.liao@mail.utoronto.ca)
% on March 15, 2016.

% calculate the drift
nu = mu - sigma * sigma/2;

% generate potential paths
% S = S0 * [ones(1, nsims); ...
%           cumprod(exp(nu * dt + sigma * sqrt(dt) * randn(nsteps, nsims)), 1)];

% transpose is what I want
S = transpose(S0 * [ones(1, nsims); ...
              exp(cumsum(nu * dt + sigma * sqrt(dt) * randn(nsteps, nsims)))]);

% ^^ I think these are interchangeable.

end

