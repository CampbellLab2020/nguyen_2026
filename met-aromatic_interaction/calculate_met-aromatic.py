## Calculate the angle and distance between methionine sulfur atom
## and center of aromatic residue in a5b1:biig2 MD simulations.
## Requires an csv file of residue atom indicies

import numpy as np
import glob
import pandas as pd
import mdtraj as md
import os

## Calculate angle from center of ring to sulfur atom
def measure_ring_sulf_angle(t, ring1_atoms, sulf_atom):
    angle = []
    
    for frame in range(t.n_frames):
        ring1_xyz = t.xyz[frame, ring1_atoms]
        ring1_xyz_center = np.sum(t.xyz[frame,ring1_atoms], axis=0)/3
        
        v1 = ring1_xyz[1] - ring1_xyz[0]
        v2 = ring1_xyz[2] - ring1_xyz[0]
        n1 = np.cross(v1,v2)
        n1 = n1/np.linalg.norm(n1) #normalize to unit vector

        sulf_xyz = t.xyz[frame, sulf_atom]
        u2 = sulf_xyz - ring1_xyz_center
        n2 = u2/np.linalg.norm(u2) #normalize to unit vector


        angle1 = np.arccos(np.dot(n1,n2))
        angle2 = np.arccos(np.dot(n1,-n2))

        angle.append(np.min([angle1,angle2]))
    
    return np.array(np.degrees(angle))

## Calculate angle from center of ring to sulfur atom
def measure_ring_sulf_dist(t, ring1_atoms, sulf_atom):
    ring1_xyz_center = np.sum(t.xyz[:,ring1_atoms], axis=1)/3
    sulf_xyz = t.xyz[:, sulf_atom]
    sulf_dist = np.sqrt(np.sum(np.square(ring1_xyz_center-sulf_xyz), axis=1))
    
    return sulf_dist


biig2_df = pd.read_csv("a5b1_biig2_met-aromatic_pairs.csv")
sulf_df = biig2_df[biig2_df.type =="ring-sulf"]
print(sulf_df)

for index, row in sulf_df.iterrows():
    print(f"Measuring interaction between {row['a5b1_resn']}{row['a5b1_resid']}-{row['fab_resid']}{row['fab_resn']}")

prmtop = "../prmtop/a5b1_biig2-strip.prmtop"
traj_list = sorted(glob.glob("../strip-traj/a5b1_biig2_clone*-strip.xtc"))

for traj in traj_list:
    t = md.load(traj, top=prmtop)

    dir = os.path.dirname(traj)
    out_file = traj.replace(dir+'/','',1)
    out_file = out_file.replace("-strip.xtc","")

    for index, row in sulf_df.iterrows():
        group1 = t.topology.select(f"resid {row['resid1_mdtraj']} and name {row['atom1']}")
        group2 = t.topology.select(f"resid {row['resid2_mdtraj']} and name {row['atom2']}")[0]

        angle = measure_ring_sulf_angle(t, ring1_atoms=group1, sulf_atom=group2)
        dist = measure_ring_sulf_dist(t, ring1_atoms=group1, sulf_atom=group2)

        np.save(f"{out_file}-{row['a5b1_resn']}{row['a5b1_resid']}-{row['fab_resid']}{row['fab_resn']}-angle.npy", angle)
        np.save(f"{out_file}-{row['a5b1_resn']}{row['a5b1_resid']}-{row['fab_resid']}{row['fab_resn']}-dist.npy", dist)

    print(traj)
