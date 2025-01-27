# SimpleFEA - A super-simple 2D FEA implementation

This repository contains a super-simple 2D finite-element analysis implementation for solving the textbook linear div-grad equation with 1st-order triangular elements.

It is inefficient, non-vectorized, and will gorge up all available RAM on problems with more than a few thousand nodes (due to the stiffness matrix being assembled as non-sparse), but it should be as easy
to understand as possible for a relatively low-level FEA implementation.

