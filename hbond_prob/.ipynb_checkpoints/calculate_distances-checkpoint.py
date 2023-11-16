import os
import numpy as np
import mdtraj as md
import pandas as pd
import glob 
import itertools

def measure_ring_ring_dist(t, ring1_atoms, ring2_atoms):
    ring1_xyz_center = np.sum(t.xyz[:,ring1_atoms], axis=1)/3
    ring2_xyz_center = np.sum(t.xyz[:,ring2_atoms], axis=1)/3

    ring_dist = np.sqrt(np.sum(np.square(ring1_xyz_center-ring2_xyz_center), axis=1))

    return ring_dist

def calculate_ring_ring_percent(t, ring1_atoms, ring2_atoms):
    dist = measure_ring_ring_dist(t, ring1_atoms, ring2_atoms)

    count = np.zeros(len(dist))
    for i in range(len(dist)):
        if dist[i] < 0.70:
            count[i] = 1

    freq = np.sum(count)/len(count)

    return freq

def calculate_min_dist_hbond_percent(t, group1, group2, cutoff = 0.35):
    atom_list = itertools.product(group1, group2)
    atom_list = np.array(list(atom_list))
    
    dist = md.compute_distances(t, atom_list)
    
    #find min distance in each frame
    min_dist = np.min(dist, axis=1) 
    
    #count instances min_dist is less the cutoff
    hbond_formed = np.sum(min_dist < cutoff)/len(min_dist)
    
    return hbond_formed


######################################
### a5b1:BIIG2 interface distances ###
######################################

parm = "../prmtop/a5b1_biig2-strip.prmtop"
system = parm.replace("-strip.prmtop","")
last_frames = 20000
first_frame = 500
biig2_df = pd.read_csv("a5b1_biig2_dist_pairs.csv")
traj_list = sorted(glob.glob("../strip-traj/a5b1_biig2_clone*-strip.xtc"))

for clone in traj_list:
    dir = os.path.dirname(clone)
    clone = clone.replace(dir+'/','',1)
    out_file = clone.replace("-strip.xtc","")

    t = md.load(f"../strip-traj/{clone}", top=parm)[first_frame:last_frames]
    freq = []

    ##Calculate hbond probabilities
    for index, row in biig2_df.iterrows():
        group1 = t.topology.select(f"resid {row['resid1_mdtraj']} and name {row['atom1']}")
        group2 = t.topology.select(f"resid {row['resid2_mdtraj']} and name {row['atom2']}")


        hbond_percent = calculate_min_dist_hbond_percent(t, group1, group2, cutoff=0.35)
        freq.append(hbond_percent)

    biig2_df[f"{out_file}_freq"] = freq

    print(f"{clone}, {t.n_frames}")

biig2_df.to_csv("a5b1_biig2_hbond_freq.csv")


#####################################
### a5b1:MINT interface distances ###
#####################################

parm = "../prmtop/a5b1_mint-strip.prmtop"
system = parm.replace("-strip.prmtop","")
last_frames = 20000
first_frame = 500
traj_list = sorted(glob.glob("../strip-traj/a5b1_mint_clone*-strip.xtc"))
mint_df = pd.read_csv("a5b1_mint_dist_pairs.csv")

for clone in traj_list:

    dir = os.path.dirname(clone)
    clone = clone.replace(dir+'/','',1)
    out_file = clone.replace("-strip.xtc","")
       
    t = md.load(f"../strip-traj/{clone}", top=parm)[first_frame:last_frames]
    
    freq = []

    ##Calculate hbond probabilities
    hbonds_df = mint_df[mint_df.type =="hbond"]
    for index, row in hbonds_df.iterrows():
        group1 = t.topology.select(f"resid {row['resid1_mdtraj']} and name {row['atom1']}")
        group2 = t.topology.select(f"resid {row['resid2_mdtraj']} and name {row['atom2']}")

        hbond_percent = calculate_min_dist_hbond_percent(t, group1, group2, cutoff=0.35)
        freq.append(hbond_percent)

    ##Calculate ring_ring probabilities
    ring_df = mint_df[mint_df.type =="ring-ring"]
    for index, row in ring_df.iterrows():
        group1 = t.topology.select(f"resid {row['resid1_mdtraj']} and name {row['atom1']}")
        group2 = t.topology.select(f"resid {row['resid2_mdtraj']} and name {row['atom2']}")

        ring_percent = calculate_ring_ring_percent(t, ring1_atoms=group1, ring2_atoms=group2)
        freq.append(ring_percent) 

    mint_df[f"{out_file}_freq"] = freq

    print(f"{clone}, {t.n_frames}")

mint_df.to_csv("a5b1_mint_hbond_freq.csv")