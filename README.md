# PhasePack
Phase retrieval is the recovery of a signal from only the magnitudes, and not the phases, of complex-valued linear measurements.  Phase retrieval problems arise in many different applications,
particularly in crystallography and microscopy.

PhasePack is a collection of sub-routines for solving classical phase retrieval problems.  PhasePack contains implementations of both classical and contemporary phase retrieval routines.  It also contains scripts that apply these signal recovery routines to datasets, and compares the results of different methods in terms of speed and quality.  PhasePack can do comparisons not only with synthetic datasets, but also with several publicly available real-world datasets.

The routines implemented in PhasePack share a simple common interface, making it easy to switch between
solvers.  The interface also gives users the option to control the initialization and runtime of different algorithms, making it easy to "plug and play" different routines and options.

# To find out more
See the main [PhasePack page](http://cs.umd.edu/~tomg/projects/phasepack/) for complete information, or check out the [user guide](https://arxiv.org/abs/1711.09777).
