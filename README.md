# MD
This repository contains a **from-scratch molecular dynamics (MD) engine written in C++** to demonstrate MD algorithms and C++ coding proficiency. It numerically integrates [Newton’s equations](https://en.wikipedia.org/wiki/Newton%27s_laws_of_motion) for particles with bonded and non-bonded interactions using neighbor lists and periodic boundary conditions. Bonded terms include harmonic bonds, angles, and dihedrals; non-bonded terms include the [Lennard-Jones 12-6 potential](https://en.wikipedia.org/wiki/Lennard-Jones_potential) (for van der Waals forces) and [Coulombic electrostatics](https://en.wikipedia.org/wiki/Coulomb%27s_law). Long-range Coulomb forces are handled by the [Particle Mesh Ewald (PME)](https://en.wikipedia.org/wiki/Ewald_summation#Particle_mesh_Ewald_(PME)_method) method, which splits interactions into real-space and Fourier-space parts for rapid convergence.

It is a personal project built to demonstrate understanding of **molecular dynamics algorithms**, **numerical methods**, and **modern C++ design**, rather than a production-ready MD package.

The code implements core components found in real MD engines (bonded and non-bonded forces, neighbor lists, PME electrostatics, FFTs, periodic boundary conditions) with an emphasis on correctness, clarity, and performance-aware design.

> This project is **not intended for real scientific use**. It is a learning and portfolio project.

---

## Features

### Force Field Components

- **Bonded interactions**
  - Harmonic bond stretching
  - Angle bending
  - Dihedral (torsional) potentials

- **Non-bonded interactions**
  - Lennard-Jones 12-6 potential
  - Coulomb electrostatics

### Long-Range Electrostatics

- **Particle Mesh Ewald (PME)**
  - Real-space short-range electrostatics
  - Reciprocal-space electrostatics using FFTs
  - Charge spreading onto a mesh and force interpolation

- **FFT-based implementation**
  - Uses FFTW3 for fast Fourier transforms
  - O(N log N) scaling for long-range electrostatics

### Performance-Oriented Design

- **Neighbor / pair lists**
  - Verlet-style pair list construction
  - Reduces cost of short-range non-bonded interactions

- **Periodic boundary conditions**
  - Minimum image convention
  - Supports fully periodic simulation boxes

---
