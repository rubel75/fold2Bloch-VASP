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
Perform the band structure calculation and generate WAVECAR file for the k-mesh of interest. Then execute

`/path/to/fold2Bloch WAVECAR FX:FY:FZ [-ncl]`

Options:

  `WAVECAR` -- name of the input VASP wavefunction file

  `FX:FY:FZ` -- multiplicity of the primitive cell that was used to construct a supercell

  `-ncl` -- the WAVECAR is produced by vasp_ncl code, which implies that the wavefunctions are spinors (defauls: assumption that WAVECAR comes from vasp_std or vasp_gam)

Examples of "real life" applications can be found in [arXiv:1405.4218](http://arxiv.org/abs/1405.4218)

Please communicate your feedback, support or feature requests via WIEN2k [mailing list](http://www.wien2k.at/reg_user/mailing_list)
