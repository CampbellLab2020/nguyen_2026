# Structural and functional characterization of integrin &alpha;5-targeting antibodies for anti-angiogenic therapy

Contains simulation files and analysis scripts for processing and generating figures related to MD simulations of &alpha;5&beta;1:Fab complex.

1. `getcontacts_interaction/` - uses the `getcontacts.io` package to identify contacts.
1. `hbond_prob/` - distance calculation between hydrogen bond acceptor and donor atoms.
1. `met-aromatic_interaction/` - distance and angle calculation between methionine sulfur atoms and aromatic residues.
1. `prmtop` - AMBER topology and parameter for production simulations and analysis.
1. `ref` - reference structure file in `xtc` format.
1. `rmsd-rmsf-CA` - RMSD and RMSF calculation using the `CPPTRAJ` module.