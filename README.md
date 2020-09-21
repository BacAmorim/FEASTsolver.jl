# FEASTsolver

A Julia implementation of the FEAST eigensolver algorithm [1]. Being develloped for learning and testing purposes.

Current functionality:
* Cauchy contour integral implemented using Gauss-Legendre quadrature along a circle.
* Solution of linear system using direct solver based on LU decomposition.

Limitation:
* The subspace dimension must estimated.

## References
[1] "A Density Matrix-based Algorithm for Solving Eigenvalue Problems", Eirv Polizzi, Phys. Rev. B 79, 115112 (2009), e-print: arXiv:0901.2665[](https://arxiv.org/abs/0901.2665)

