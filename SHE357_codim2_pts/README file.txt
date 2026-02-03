README file

Matlab Scripts for numerical computations and continuation of codimension two points for the cubic-quintic-septic Swift-Hohenberg equation.

This library can do the following: 

1. Compute families of periodic orbits and the corresponding spatial Hamiltonian. To do this, run: run_SHE357_periodic.m 

2. Continue a branch of periodic solutions where the spatial Hamiltonian possesses a critical point. This is done by running run_SHE357_Hk_0.m
The code assumes that a branch of periodic solutions has already been computed, and a critical point of the spatial Hamiltonian has been identified. 

3. Continue branch of periodic solutions at the codimension two point such that H_k = H_kk = 0. To do this run run_SH357_codim2_point.m 
A good initial guess from the family of solutions computed in code 2 is required for reasonably fast convergence. 

Continuation uses code developed by Danielle Avitabile (SecantContinuation.m). 
A small modification was made since the spectrum of the periodic orbits correspond to its Floquet multipliers. 