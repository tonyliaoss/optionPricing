# Option Pricing

A comparison of different numerical methods used for option pricing, for both
the pricing of European and American options. We are interested in the put
options because the prices between European and American put options differ,
whereas American call options always have the same value as European call
options through a simple arbitrage argument.

Basically we want to solve the Black-Scholes PDE, which is a second-order,
linear, backward-parabolic, partial differential equation. I have done some
parametrization and change of variables to make the math slightly easier, which
I will explain later...

