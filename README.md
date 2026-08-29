# Shared ligand-blocking mechanism but distinct conformational modulation by α5-targeting antibodies BIIG2 and MINT1526A

**Adam Nguyen<sup>1,2</sup>, Joel B. Heim<sup>3,4</sup>, Gabriele Cordara<sup>4</sup>, Matthew C. Chan<sup>1</sup>, Hedda Johannesen<sup>4</sup>, Cristine Charlesworth<sup>5</sup>, Ming Li<sup>3</sup>, Caleigh M. Azumaya<sup>1,6</sup>, Benjamin Madden<sup>5</sup>, Ute Krengel<sup>4</sup>, Alexander Meves<sup>3</sup>, Melody G. Campbell<sup>1,2,7</sup>**

<sup>1</sup> Basic Sciences Division, Fred Hutchinson Cancer Center, Seattle, WA, 98109, United States of America

<sup>2</sup> Biological Physics Structure and Design Program, University of Washington, Seattle, WA, 98195, United States of America

<sup>3</sup> Department of Dermatology, Mayo Clinic, Rochester, MN, 55905, United States of America

<sup>4</sup> Department of Chemistry, University of Oslo, Oslo, 0315, Norway

<sup>5</sup> Medical Genome Facility, Proteomics Core, Mayo Clinic, Rochester, MN, 55905, United States of America

<sup>6</sup> Present address: Genentech, South San Francisco, CA 94080, USA

<sup>7</sup> Lead contact

*Structure* 2026. DOI: [https://10.1016/j.str.2026.08.001](https://doi.org/10.1016/j.str.2026.08.001)

---
This repository contains simulation files and analysis scripts for processing and generating figures related to MD simulations of &alpha;5&beta;1:Fab complex.

1. `getcontacts_interaction/` - uses the `getcontacts.io` package to identify contacts.
1. `hbond_prob/` - distance calculation between hydrogen bond acceptor and donor atoms.
1. `met-aromatic_interaction/` - distance and angle calculation between methionine sulfur atoms and aromatic residues.
1. `prmtop/` - AMBER topology and parameter for production simulations and analysis.
1. `ref/` - reference structure file in `xtc` format.
1. `rmsd-rmsf-CA/` - RMSD and RMSF calculation using the `CPPTRAJ/` module.
1. `a5b1_antibody_MD_analysis.ipynb` - python notebook containing scripts to generate associated figures.
