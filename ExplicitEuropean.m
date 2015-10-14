function [ P ] = ExplicitEuropean( S, tau, E, r, sigma )
% ExplicitEuropean: Computes the fair value of the European put option
% using the explicit Euler method.
%
% Notes
%   - This algorithm does not solve the Black-Scholes PDE explicitly; it
%   tries to solve a parametrized version of the Black-Scholes PDE. I have
%   made the change of variables tau = T - t (which converts the backwards
%   parabolic equation to a forward parabolic one), and S = S_0 * exp(x),
%   (which gets rid of a bunch of constants in front of the d^2/dS^2 term)
%   So I'm solving the equation:
%      dP/dtau = sigma^2/2 * d^2 P/ dx^2 + r * dP/dx - rP
%   - alpha must be kept below ~0.5 for stability.
%
% INPUT PARAMETERS
%   S - the current market value of the asset.
%   tau - time to expiry of the option, expressed in years.
%   E - the strike price, or exercise price, of the option.
%   r - the risk-free investment return per annum.
%   sigma - the market volatility.

alpha = 0.4; % alpha = dt/(dx^2)
numPartitionsX = 500;

% determine the maximum and minimum x values.
z = 3; % variable for random walk
S_t = S * exp((r - 0.5 * sigma^2) * tau + sigma * sqrt(tau) * z);
x_max = log(S_t / S);
x_min = -x_max;

dx = (x_max - x_min) / numPartitionsX;
dt = alpha * dx^2;
numPartitionsT = ceil(tau/dt);

% we first populate PRICE, a row vector
PRICE = zeros(numPartitionsT + 1, numPartitionsX + 1);
% then create X, a row vector
X = linspace(x_min, x_max, numPartitionsX + 1);
% create T, a row vector
T = linspace(0, tau, numPartitionsT + 1);

% initial conditions
PRICE(1,:) = max((E - S * exp(X)), 0);

% boundary conditions
% as x approaches -inf
PRICE(:, 1) = E * exp(-r * T); % present value of payoff
% as x approaches inf
PRICE(:, end) = 0;

% populate the payoff matrix
for t = 2:numPartitionsT+1
    for x = 2 : numPartitionsX
        PRICE(t, x) = sigma^2 * alpha / 2 * (PRICE(t-1, x+1) - 2 * PRICE(t-1, x) + PRICE(t-1, x-1)) ...
                       + r * dt / dx * (PRICE(t-1, x+1) - PRICE(t-1, x)) ...
                       + (1 - r * dt) * PRICE(t-1, x);
        if isnan(PRICE(t,x)) % error: method diverged
            % catch error
            fprintf('Diverged at: t = %i, x = %i\n', t, x);
            fprintf('%f %f %f %f\n', PRICE(t-1, x+1), PRICE(t-1, x), PRICE(t-1, x-1), PRICE(t,x));
            return
        end
    end
end

% return the 'middle' value of X
P = PRICE(numPartitionsT+1, ceil((numPartitionsX + 1) / 2));

end
