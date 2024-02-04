Better check for contractions in newton raphson method: |(xk - r)| < |(xk+1 - r)| for example
Bigger Test Suite
Imlement backward error check:
    1. backward error (relative and absolute): Slides Lecture 9 page 16
    2.
Implement condition number calculation
Implemet different accuracy at different stages of calculation based on condition number
Condition number of A := kp(A)
E.g.: kp(A) = 10^n then n digits will be lost in solving Ax = b
1 < kp(A) < 10^8 - well conditioned
10^8 <= kp(A) < 10^16 - ill conditioned
10^16 <= kp(A)  - completely unreliable
kp(a) = ||A||p||A^-1||p - expensive to calculate

Implement Gauss and Lu factorization by pivoting (use the custom pointer initilization matrix in Jonas Skeppstedts book )
