function [Z,E] = solve_lrr(X,lambda)
Q = orth(X');
A = X*Q;
[Z,E] = lrra(X,A,lambda,false);
Z = Q*Z;