# Tree-Tensor-Network States for Bosons
This repository contains a Density Matrix Renormalization Group (DMRG) implementation based on Matrix Product States (MPS), which is capable of simulating arbitrary bosonic tree-structured Hamiltonians. It implements the Optimized Boson Basis (OBB) and a variational shift to allow for almost arbitrary local bosonic Hilbert space dimensions, allowing to access the strong coupling regime.
The time-evolution is performed with the Time-Dependent Variational Principle (TDVP), based on the Dirac-Frenkel Principle.

If you intend to use this code, please cite the following publications which this work is based upon:

1. C. Guo, A. Weichselbaum, J. von Delft, and M. Vojta, Phys. Rev. Lett. 108, 160401 (2012). DOI: [10.1103/PhysRevLett.108.160401][Guo2012]  
2. F. A. Y. N. Schröder and A. W. Chin, Phys. Rev. B 93, 75105 (2016). DOI: [10.1103/PhysRevB.93.075105][Schroeder2016]  


## Key features
- 1D Hamiltonians, e.g. Spin-Boson Model
- Hamiltonians with **Tree topology**
- Bosons with local Hilbert space **~ 200**
- Fully **dynamic** bond dimensions
- Ground state calculations
- Time-evolution with **TDVP**
- Supports Open Quantum Systems with **arbitrary bosonic environments**
- Extraction of **Memory Kernel**/Transfer Tensors (weak coupling)
- Transfer Tensor Method (**TTM**) for rapid extrapolation to **long times**
- Adiabatic Total Energy Surfaces (**TES**) and Potential Energy Surfaces (**PES**)


## Publications using this implementation:
[1] F. A. Y. N. Schröder and A. W. Chin, Phys. Rev. B 93, 75105 (2016). DOI: [10.1103/PhysRevB.93.075105][Schroeder2016]  

## This program is based on:
|Technique|Paper|
|---|---|
|**OBB**|C. Guo, A. Weichselbaum, J. von Delft, and M. Vojta, Phys. Rev. Lett. 108, 160401 (2012). DOI: [10.1103/PhysRevLett.108.160401][Guo2012]|
| |C. Zhang, E. Jeckelmann, and S. R. White, Phys. Rev. Lett. 80, 2661 (1998). DOI: [10.1103/PhysRevLett.80.2661][Zhang1998] |
|**TDVP**|J. Haegeman, C. Lubich, I. Oseledets, B. Vandereycken, and F. Verstraete, Phys. Rev. B 94, 165116 (2016). DOI: [10.1103/PhysRevB.94.165116][Haegeman2016]|
|**EXPOKIT**|R. B. Sidje, ACM Trans. Math. Softw. 24, 130 (1998). DOI: [10.1145/285861.285868][Sidje1998]|
|**TTM**|J. Cerrillo and J. Cao, Phys. Rev. Lett. 112, 1 (2014). DOI: [10.1103/PhysRevLett.112.110401][Cerillo2014]|

[Guo2012]: https://doi.org/10.1103/PhysRevLett.108.160401
[Schroeder2016]: https://doi.org/10.1103/PhysRevB.93.075105
[Sidje1998]: https://doi.org/10.1145/285861.285868
[Zhang1998]: https://doi.org/10.1103/PhysRevLett.80.2661
[Haegeman2016]: https://doi.org/10.1103/PhysRevB.94.165116
[Cerillo2014]: https://doi.org/10.1103/PhysRevLett.112.110401
