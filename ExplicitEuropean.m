function [ P ] = ExplicitEuropean( S, tau, E, r, sigma )
% Computes the fair value of the European put option using the explicit
% finite differences method.
%   Detailed explanation goes here

alpha = 0.2; % alpha = dt/(dx^2)
numPartitionsX = 100;

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
T = linspace(0, tau, numPartitionsT + 1);

% initial conditions
for i = 1:numPartitionsX
    PRICE(1,i) = max((E - S * exp(X(i))), 0);
end

% boundary conditions
for i = 1:numPartitionsT+1
    % as x approaches -inf
    PRICE(i,1) = E * exp(-r * T(i));
    % as x approaches inf
    PRICE(i, numPartitionsX+1) = 0;
end

% populate the payoff matrix
for t = 2:numPartitionsT+1
    for x = 2 : numPartitionsX
        PRICE(t, x) = sigma^2 * alpha / 2 * (PRICE(t-1, x+1) - 2 * PRICE(t-1, x) + PRICE(t-1, x-1)) ...
                       + r * dt / dx * (PRICE(t-1, x+1) - PRICE(t-1, x)) ...
                       + (1 - r * dt) * PRICE(t-1, x);
        if isnan(PRICE(t,x))
            fprintf('Diverged at: t = %i, x = %i\n', t, x);
            fprintf('%f %f %f %f\n', PRICE(t-1, x+1), PRICE(t-1, x), PRICE(t-1, x-1), PRICE(t,x));
            return
        end
    end
end

P = PRICE(numPartitionsT+1, ceil((numPartitionsX + 1) / 2));

end
