# Cerisse

**A research computational solver for Large Eddy Simulation  and Direct Numerical Simulation  of high-speed compressible turbulent reactive flows.**

Cerisse is a high-performance CFD framework for research in compressible reacting flows, designed for modern heterogeneous HPC architectures using MPI and GPU acceleration.

## Capabilities

- High-order finite-volume solver
- Large Eddy Simulation (LES) and Direct Numerical Simulation (DNS)
- Parallel execution with MPI + GPU acceleration via AMReX
- Adaptive Mesh Refinement (AMR)
- Complex geometries via Immersed Boundaries (IBM) or Embedded Boundaries (EBM)
- Flexible thermodynamics and detailed chemistry integration through PelePhysics
- Modular research-oriented architecture for custom physics development

## Quick Install

Before installing, ensure the required dependencies are available (C++ compiler, MPI, AMReX, optional PelePhysics support, and GPU toolchains as needed).

```bash
cd cerisse/lib
./install.sh safe
```

For full dependency requirements and installation options, see the documentation website (or run script `bin/checkreq.sh`)

## Documentation

Full installation instructions, requirements, tutorials, and user documentation are available at:

**https://hslesdnsrf.gitbook.io/cerisse-docs**

## Getting Started

1. Check requirements (documentation/script)
2. Install dependencies
3. Run the quick install command above
4. Compile and run the tutorial cases in `tst/tutorial`

## License

BSD 3-Clause License
