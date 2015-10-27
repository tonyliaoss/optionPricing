The system of linear equations that we are solving is the same as in "simple",
but instead of using Matlab's built-in `\` operator to solve the system of
linear equations, I define a custom method based on the **successive
over-relaxation** (SOR) technique.

The SOR is an iterative technique for solving systems of linear equations, and
this makes them easily extendable to computing the prices of American options,
which have no simple, closed-form solutions.

