#!/usr/bin/env python


from __future__ import print_function
import mdtraj as md
import numpy as np
from scipy.spatial.distance import cdist, squareform


def pdb_to_nmr(resid):
    '''
    Convert from simulated PDB to NMR indexing.

    The resulting index will match what is in `trajectory.pdb`. We use
    the convention that a residue not present in the PDB file is `None`.

    Restype		NMR Index		PDB Index		Note
    ALA         1               None            First protein residue
    ASP         2               None
    GLN         3               None
    LEU         4               None
    THR         5               1
    ...
    THR         146             142
    ALA         147             None
    LYS	        148             None
    CYS         149             None            Last protein residue

    CA
    CA
    CA
    CA

    ALA         201             None            First peptide residue
    ARG         202             147
    ARG         203             148
    ...
    LEU         218             163
    SER         219             164
    SER         220             165             Last peptide residue
    '''
    resid = int(resid)
    if resid >= 1 and resid <= 142:
        return resid + 4
    elif resid >= 147 and resid <= 165:
        return resid + 55
    else:
        raise RuntimeError()


traj = md.load('trajectory.pdb')

nitrogens = traj.topology.select('name N')
spin_labels = traj.topology.select('name OND')

distances = np.zeros((len(nitrogens), len(spin_labels), traj.n_frames))
for i in range(traj.n_frames):
    distances[:, :, i] = cdist(traj.xyz[i, nitrogens, :], traj.xyz[i, spin_labels, :])

# do 1/r^6 average
avg_distances = np.mean(distances**(-6.), axis=2)**(-1.0 / 6.0)


with open('average_distances.dat', 'w') as outfile:
    for j, spin_atom_ind in enumerate(spin_labels):
        for i, N_atom_ind in enumerate(nitrogens):
            # +1 is to convert from zero-based indexing of mdtraj
            # to 1-based index of PDB file
            N_res = pdb_to_nmr(traj.topology.atom(N_atom_ind).residue.index + 1)
            res_type = traj.topology.atom(N_atom_ind).residue.name
            spin_res = pdb_to_nmr(traj.topology.atom(spin_atom_ind).residue.index + 1)

            print('{}\t{}\t{}\t{}'.format(spin_res, N_res, res_type,
                                          avg_distances[i, j]), file=outfile)
