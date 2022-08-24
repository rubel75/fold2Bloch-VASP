## fold2Bloch

The fold2Bloch utility is designed to unfold the band structure of a supercell obtained with the Vienna Ab initio Simulation Package ([VASP](https://www.vasp.at)) and compute an effective band structure in a primitive representation. It facilitates interpretation of large-scale electronic structure calculations, where the Bloch character of electronic eigenstates is perturbed by a disorder (defects, alloy elements, etc). Reading of wavefunctions is adapted from the [WaveTrans](http://www.andrew.cmu.edu/user/feenstra/wavetrans) code.

### Contributors:
* Oleg Rubel (McMaster University <rubelo@mcmaster.ca>)
* Michael Widom and Randall Feenstra (Carnegie Mellon University)
* Anton Bokhanchuk (Confederation College)

### Installation:
First clone the GitHub repository

`$ git clone https://github.com/rubel75/fold2Bloch-VASP`

The `makefile` is set up for Intel Fortran compiler `ifort`. To compile, simply execute

`$ cd fold2Bloch-VASP; make`

### Execution
First, you need to perform the band structure calculation and generate WAVECAR file for the k-mesh of interest using standard VASP procedure (see VASP [guidelines](https://cms.mpi.univie.ac.at/wiki/index.php/Si_bandstructure)). (The following Matlab script `utils/fold.m` is designed to assist with preparing a folded string of k-points that will unfold on a desired k-path.) It is advised to increase the number of empty bands (NBANDS=... in INCAR file) by a factor of 1.2-2 beyond a VASP-proposed default value to get a resonable description of higher energy states.

Once WAVECAR is ready, execute

`/path/to/fold2Bloch WAVECAR "P11 P12 P13:P21 P22 P23:P31 P32 P33" [-ncl]`

Options:

  `WAVECAR` -- name of the input VASP wavefunction file

  `"P11 P12 P13:P21 P22 P23:P31 P32 P33"` -- transformation matrix from primitive lattive verctors `a_p` to supercell lattive verctors `a_s` (same as in [VESTA](https://jp-minerals.org/vesta/en/) or the [Bilbao Crystallographic Server](https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-doc-trmat)):

  ```
  a_s(i) = sum_j a_p(j)*P(j,i)      i,j = 1, 2, 3
  ```

  `-ncl` -- optional switch that needs to be activated when the WAVECAR is produced by vasp_ncl code, which implies that the wavefunctions are spinors (default assumption is that WAVECAR comes from vasp_std or vasp_gam).

### Output
Output is written to `WAVECAR_*.f2b` file(s). There is one output file for non-spin-polarized calculation and two files for the spin-polarized or spinor (vasp_ncl) calculations. Below is a sample of an output file.

    New K-values (x, y, z)       Eigenvalue (eV)  Weight
    0.000000   0.000000   0.000000 -22.019998   1.000000
    0.000000   0.200000   0.000000 -22.019998   0.000000
    0.000000   0.400000   0.000000 -22.019998   0.000000
    0.000000  -0.400000   0.000000 -22.019998   0.000000
    0.000000  -0.200000   0.000000 -22.019998   0.000000
    0.200000   0.000000   0.000000 -22.019998   0.000000
    ...

### Plotting results
The Matlab code `utils/ubs_dots_VASP.m` or the octave code `utils/ubs_bmp.m` are designed to assist with plotting the band structure and the Bloch weights. Please refer to its input section for description of user input variables. Each plotting tool has its own pros/cons. Below is a sample of the Matlab plot for the unfolded band structure of a dilute GaP:N alloy. It is followed by a different style plot of a dynamic band structure in a perovskite lattice generated using Octave `utils/ubs_bmp.m` and gnuplot `utils/f2b-band-structure.plt` from [Phys. Rev. Materials **2**, 114604 (2018)](https://doi.org/10.1103/PhysRevMaterials.2.114604).

<img src="https://github.com/rubel75/fold2Bloch-VASP/blob/master/graphics/GaP%2BN.png" width="600">

<img src="https://github.com/rubel75/fold2Bloch-VASP/blob/master/graphics/perovskite.png" width="350">


### References
If you find the results useful and publishable, we will appreciate citing the following papers:
* O. Rubel, A. Bokhanchuk, S. J. Ahmed, and E. Assmann "Unfolding the band structure of disordered solids: from bound states to high-mobility Kane fermions", [Phys. Rev. B **90**, 115202 (2014)](http://olegrubel.mcmaster.ca/publications/2014/Rubel_PRB_90_115202.pdf).
* L.-W. Wang, L. Bellaiche, S.-H. Wei, and A. Zunger "Majority representation of alloy electronic states", [Phys. Rev. Lett. **80**, 4725 (1998)](https://doi.org/10.1103/PhysRevLett.80.4725).
