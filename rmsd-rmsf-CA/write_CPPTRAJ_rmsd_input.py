## This script writes cpptraj input scripts to calculate
## C-alpha RMSD and RMSF for the simulations.
## For both calculations, the initial structure for 
## simulation is used as the reference structure.

import glob
import os


for parm in glob.glob("../prmtop/a5b1*prmtop"):
    dir = os.path.dirname(parm)
    system = parm.replace(dir+'/','',1)
    system = system.replace("-strip.prmtop","")
    ref = f"{system}-tleap-ref.xtc"
    
    for clone in glob.glob(f"../strip-traj/{system}_clone*xtc"):
        
        dir = os.path.dirname(clone)
        clone = clone.replace(dir+'/','',1)
        out_file = clone.replace("-strip.xtc","")
        file = open(f"cpp_{clone}-analysis.in","w")
        file.write(f"parm {parm} \n")
        file.write(f"trajin ../strip-traj/{clone} \n")
        file.write(f"reference ../ref/{ref} [ref] \n")
        file.write(f"rms @CA reference [ref] \n")
        file.write(f"atomicfluct out {out_file}-rmsf-CA.dat :1-697@CA byres \n")
        file.write(f"rmsd @CA ref [ref] out {out_file}-rmsd-CA.dat \n")
        file.write("run")
