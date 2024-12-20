## DL_POLY-Quantum-v2.1

DL_POLY Quantum is a general purpose open-source software package for classical and advanced path integral simulations in condensed phases subject to nuclear quantum effects.
Developed by:
    Nathan London 
    Dil K. Limbu 
    Md Omar Faruque
    Mohammad R. Momeni
    and Farnaz A. Shakib 

The original DL_POLY Quantum, version 1.0, which was developed by Mohammad R. Momeni and Farnaz A. Shakib at the New Jersey Institute of Technology. DL_POLY Quantum v1.0, in turn, was a modified version of DL_POLY Classic v1.10 originally developed in the Daresbury Laboratory.

The documentaion website for DL_POLY Quantum can be found at https://dlpolydocs.readthedocs.io . The documentation includes details for compiling the code as well as tutorials for performing different types of simulations available in the software

For questions and feedbacks please contact Dr. Momeni (mmomenitaheri@umkc.edu) or Dr. Shakib (shakib@njit.edu).

## Version 2.1 Updates

- Inclusion of support for f-QCMD simulations through the addition of the fqcmd_module
- Inclusion of support for f-CMD using iterative Boltzmann inversion (similar to f-QCMD)
- Inclusion of support for h-CMD simulations for complex, heterogeneous systems
- Update to the correlation_module to calculate the total cell dipole moment instead of molecular dipoles
- Update to spectra caclulation to use dipole derivative and block average over trajectories
- Update to initial momentum distribution of PI simulations to improve accuracy of momentum distribution
- Updated PI simualations to remove center of mass drift of the system to improve accuracy of certain dynamical properties like diffusion coefficients
- Updated massive Nose-Hoover chain thermostat to work for PIMD simulations with a single bead
