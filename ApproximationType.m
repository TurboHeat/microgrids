classdef ApproximationType < uint8
%{
0: No approximation - might be slow.
1: Additive approximation - solution costs no more than epsilon + cost of optimal solution.
2: Multiplicative approximation - solution costs no more than (1+epsilon) * cost of the optimal solution.
3: Maximum Iterations - Maximum number of shortest-path runs done by the algorithm.
   Runs the additive approximation scheme with a corresponding epsilon.
%}
  enumeration
    Exact           (0)
    Additive        (1)
    Multiplicative  (2)
    MaxIter         (3)
  end
end
