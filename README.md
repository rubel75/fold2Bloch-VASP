## fold2Bloch

Unfolding of first-principle electronic band structure obtained with [VASP](https://www.vasp.at) DFT code. Reading of wavefunctions is adapted from the [WaveTrans](http://www.andrew.cmu.edu/user/feenstra/wavetrans) code.

### Contributors:
* Oleg Rubel (McMaster University <oleg.v.rubel@gmail.com>)
* Michael Widom (Carnegie Mellon University)
* Randall Feenstra (Carnegie Mellon University)
* Anton Bokhanchuk (Confederation College)

### Installation:
The `makefile` is set up for Intel Fortran compiler `ifort`. To compile, simply run

### Execution
Perform the band structure calculation and generate WAVECAR file for the k-mesh of interest using standard VASP procedure (see VASP [guidelines](https://cms.mpi.univie.ac.at/wiki/index.php/Si_bandstructure)). Then execute

`/path/to/fold2Bloch WAVECAR FX:FY:FZ [-ncl]`

Options:

  `WAVECAR` -- name of the input VASP wavefunction file

  `FX:FY:FZ` -- multiplicity of the primitive cell that was used to construct a supercell. Note: the supercell should be constracted on the basis of the primitive cell (not conventional).

  `-ncl` -- optional switch that needs to be activated when the WAVECAR is produced by vasp_ncl code, which implies that the wavefunctions are spinors (defauls: assumption that WAVECAR comes from vasp_std or vasp_gam)

### Output
Output is writen to `WAVECAR_*.f2b` file(s). There is one output file for non-spin-polarized calculation and two files for the spin-polarized or spinor (vasp_ncl) calculations. Below is a sample of an output file.

    New K-values (x, y, z)       Eigenvalue (eV)  Weight
    0.000000   0.000000   0.000000 -22.019998   1.000000
    0.000000   0.200000   0.000000 -22.019998   0.000000
    0.000000   0.400000   0.000000 -22.019998   0.000000
    0.000000  -0.400000   0.000000 -22.019998   0.000000
    0.000000  -0.200000   0.000000 -22.019998   0.000000
    0.200000   0.000000   0.000000 -22.019998   0.000000

### References
If you find the results usefull and publishable, we will appreciate citing the following papers:
* O. Rubel, A. Bokhanchuk, S. J. Ahmed, and E. Assmann "Unfolding the band structure of disordered solids: from bound states to high-mobility Kane fermions", [Phys. Rev. B **90**, 115202 (2014)](http://olegrubel.mcmaster.ca/publications/2014/Rubel_PRB_90_115202.pdf).
* V. Popescu and A. Zunger "Effective Band Structure of Random Alloys", [Phys. Rev. Lett. **104**, 236403 (2010)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.236403).
