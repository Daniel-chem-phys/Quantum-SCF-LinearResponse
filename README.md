# Quantum-SCF-LinearResponse
Quantum SCF & Linear Response Solver ⚛️
This repository contains the core of my work in translating high-level quantum chemistry theory into functional, high-performance code. Written in Fortran, this project started as a basic Self-Consistent Field (SCF) engine and eventually grew into a specialized tool for analyzing how molecules respond to external electric fields.

What's Inside
The Hartree-Fock Engine (hf.f90): This is the project's backbone. It solves the Roothaan-Hall equations through an iterative SCF loop. To keep the calculations stable and fast, I implemented Damping and DIIS (Direct Inversion in the Iterative Subspace) algorithms—essential tools for "taming" convergence in complex systems.

Linear Response & Polarizability (pol.f90): Moving beyond the ground state, I implemented a module to calculate molecular polarizability by solving Casida's equations. This part of the code handles linear response theory, shifting from simple energy calculations to full polarizability tensors.

Technical "Battle Scars"
HPC & Linear Algebra: The code is built on top of LAPACK and BLAS. You’ll see heavy use of routines like DGESV for linear systems, DGEMM for matrix products, and DSYEV for diagonalizations.

Memory Management: Everything is handled via dynamic allocation, ensuring efficient memory use for large-scale tensors and matrices.

A Note on the Environment: These files represent the algorithmic "brain" I developed within a high-performance academic infrastructure. While the support modules (utils.f90 and hfmod.f90) are not included due to university access restrictions, the core logic and numerical implementation are fully transparent and documented here.

Why This Code?
I’m sharing this because it reflects my ability to take abstract theory and turn it into robust, fast-running algorithms. It’s a mix of scientific rigor, heavy-duty linear algebra, and numerical optimization—the kind of work I truly enjoy.
