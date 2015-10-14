function  [ P ] = ImplicitEuropean( S, tau, E, r, sigma )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
numPartitionsX = 500;
numPartitionsT = 500;

r_days = r / 365;
z = 3; % variable for random walk
S_t = S * exp((r_days - 0.5 * sigma^2) * tau + sigma * sqrt(tau) * z);
x_max = log(S_t / S);
x_min = -x_max;
dx = (x_max - x_min) / numPartitionsX;
dt = tau / numPartitionsT;

PRICE = zeros(numPartitionsT + 1, numPartitionsX + 1);
X = linspace(x_min, x_max, numPartitionsX + 1);
T = linspace(0, tau, numPartitionsT + 1);

% initial condition
PRICE(1,:) = max(E - S * exp(X), 0);
% boundary condition at x = -inf
PRICE(:,1) = E * exp(-r_days * T);
% boundary condition at x = inf
PRICE(:,numPartitionsX + 1) = 0;

% coefficients
A = 0.5 * dt * sigma ^ 2 / dx^2 - 0.5 * r_days * dt / dx;
B = 1 - sigma ^ 2 * dt / dx^2 - r_days * dt/ 2 / dx;
C = 0.5 * dt * sigma^2 / dx^2 + 0.5 * r_days * dt / dx;

% generate tri-diagonal matrix
Acoeffs = zeros(numPartitionsX + 1, 1) + A;
Bcoeffs = zeros(numPartitionsX + 1, 1) + B;
Ccoeffs = zeros(numPartitionsX + 1, 1) + C;
tri = diag(Acoeffs(2:end), -1) + diag(Bcoeffs) + diag(Ccoeffs(1:end-1), +1);
tri_inv = inv(tri);

% implicit iteration
for t = 1:numPartitionsT
    temp = zeros(numPartitionsT + 1, 1);
    temp(1) = A * PRICE(t+1, 1);
    temp(end) = C * PRICE(t+1, numPartitionsT);
    b = PRICE(t,:)' -temp;
    temp = tri_inv * b;
    PRICE(t+1, (2:end-1)) = temp(2:end-1);
end

P = PRICE(numPartitionsT+1, ceil((numPartitionsX + 1)/2));

end

