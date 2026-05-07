# Convective Drying Simulation — 1D

Fortran 90 code for simulating moisture transport in porous media during 
convective drying, developed as part of the CSII course at ENSGTI/UPPA.

## Features
- 1D Finite Element Method (P1) with consistent mass matrix
- Method I: semi-implicit scheme (gravity treated implicitly)
- Method II: fractional step scheme with Runge-Kutta 2 for the hyperbolic part
- Automatic calculation of wet bulb temperature (bisection method)
- Neumann and Dirichlet boundary conditions
- Configurable time storage interval
