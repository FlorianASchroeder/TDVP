# Tree-Tensor-Network States for Bosons
This repository contains a Density Matrix Renormalization Group (DMRG) implementation based on Matrix Product States (MPS), which is capable of simulating arbitrary bosonic tree-structured Hamiltonians. It implements the Optimized Boson Basis (OBB) and a variational shift to allow for almost arbitrary local bosonic Hilbert space dimensions, allowing to access the strong coupling regime.
The time-evolution is performed with the Time-Dependent Variational Principle (TDVP), based on the Dirac-Frenkel Principle.

If you intend to use this code, please cite the following pupblications which this work is based upon:

1. C. Guo, A. Weichselbaum, J. von Delft, and M. Vojta, Phys. Rev. Lett. 108, 160401 (2012). DOI: [10.1103/PhysRevLett.108.160401][Guo2012]  
2. F. A. Y. N. Schröder and A. W. Chin, Phys. Rev. B 93, 75105 (2016). DOI: [10.1103/PhysRevB.93.075105][Schroeder2016]  


## Key features
- 1D Hamiltonians, e.g. Spin-Boson Model
- Hamiltonians with **Tree topology**
- Bosons with local Hilbert space **~ 200**
- Fully **dynamic** bond dimensions
- Ground state calculations
- Time-evolution with **TDVP**
- Compuation of Memory Kernel/Transfer Tensors
- Transfer Tensor Method (TTM) for rapid extrapolation to **long times**


## Publications using this implementation:
[1] F. A. Y. N. Schröder and A. W. Chin, Phys. Rev. B 93, 75105 (2016). DOI: [10.1103/PhysRevB.93.075105][Schroeder2016]  


[Guo2012]: https://doi.org/10.1103/PhysRevLett.108.160401
[Schroeder2016]: https://doi.org/10.1103/PhysRevB.93.075105
