# Frobenius

In the paper [provide reference], we introduce the following primality test:

* (1) Test whether $n$ is a square. If it is, declare $n$ to be composite and stop.
* (2) Set $a=1, 3, 5, \ldots$ such that the Jacobi symbol of $a^2-4$ over $n$ is $-1$, and then set $T=1, 2, 3, \ldots$ and $Q=T^2+aT+1$ such that $Q \neq |a^2-4|$ and the Jacobi symbol of $Q$ over $n$ is $-1$. If $\gcd\left((a^2-4)(a+2T)Q,n\right) \neq 1$ throughout the search, declare $n$ to be composite and stop.
* (3) If $Q^{(n-1)/2} \neq -1 \pmod{n}$, declare $n$ to be composite and stop.
* (4) If $s_n\neq -1 \pmod{n}$ or $t_n\neq a+T\pmod{n}$, declare $n$ to be composite and stop.
* (5) If $n$ is not declared composite in steps (1) to (4), declare $n$ to be a probable prime.

In this repository, we provide the underlying code and results.