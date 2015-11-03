# Projected SOR
This is a generalization of SOR for solving linear complementarity problems,
which characterize American options. Here we write out our PDE in *linear
complementarity formulation*, and then we can slightly modify our SOR algorithm
to price American options. Error checking is an interesting issue since the
pricing of American options don't have closed-form solutions.

I've checked my values against a textbook and it seems correct. An interesting
thing that you need to modify as you go from SOR to PSOR is the boundary value
of the problem: you do not discount your payoff by the present value, because
you can exercise American options at any time. That was an interesting thing
that tripped me up for a while, giving me underestimated prices for American
options.

