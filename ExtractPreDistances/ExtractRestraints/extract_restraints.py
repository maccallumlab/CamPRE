#!/usr/bin/env python


from __future__ import print_function
import pandas as pd
import numpy as np
import scipy.optimize as opt
import math
import re
import scipy.spatial.distance as sd
import random
from bokeh.plotting import figure, output_file, show
from bokeh.models.sources import ColumnDataSource
from bokeh.models import  HoverTool, Label, Span
from bokeh.layouts import gridplot


def load_r2(excel):
    df = pd.read_excel(excel, 'C149', header=None, converters={7: int})
    r2 = df.drop([0, 1, 2, 3, 4, 5, 6, 9], axis=1)
    r2.columns = ['resid', 'r2']
    r2 = r2.dropna()
    return r2


def load_data(excel, sheetname):
    df = pd.read_excel(excel, sheetname, header=None, converters={0: int})
    df.columns = [
        'resid',
        'diamagnetic',
        'c',
        'paramagnetic',
        'e',
        'f',
        'g',
        'h',
        'i',
        'j']
    # now drop the columns we don't need
    df = df.drop(['c', 'e', 'f', 'h', 'i', 'j', 'g'], axis=1)
    # now drop anything left with NaN
    df = df.dropna()
    return df


def compute_distances(para, dia, r2):
    # compute the noise level as a percentage of the highest intensity

    means = []
    stds = []
    for p, d, r in zip(para, dia, r2):
        m = compute_single_distance(p, d, r)
        means.append(m)

    return means


def compute_single_distance(para, dia, r2):
    '''
    Compute the distance from peak intensities.

    Formulas and approach from:
    W.D. Van Horn, A.J. Beel, C. Kang, C.R. Sanders, The impact of window
    functions on NMR-based paramagnetic relaxation enhancement
    measurements in membrane proteins, BBA, 2010, 1798: 140-149
    '''
    K = 1.23e-32 * 1e-12 # cm^6 s^-2 * m^6/cm^6 = m^6 s^-2
    tc = 8e-9 # ns
    omega_H = 700e6 # s^-1
    f = 4 * tc + 3 * tc / (1 + (omega_H * tc)**2 )

    ratio = para/ dia

    # Clamp ratio to lie between 1e-6 and 0.99
    ratio = 1e-6 if ratio < 1e-6 else ratio
    ratio = 0.99 if ratio > 0.99 else ratio

    func = lambda gamma, r=ratio: r2 * math.exp(-gamma * 10e-3) / (r2 + gamma) - r
    try:
        gamma = opt.newton(func=func, x0=0.1, maxiter=500)
        r = (K * f / gamma)**(1.0 / 6.0) * 1e9
        r = 5 if r > 5 else r
    except (RuntimeError, OverflowError):
        r = np.NaN
    return r


def load_crystal_distances():
    AVG_FILE = '../AverageDistances/average_distances.dat'
    df = pd.read_csv(AVG_FILE, sep='\t')
    df.columns = ['spin_label_index', 'resid', 'restype', 'distance']
    groups = df.groupby('spin_label_index')
    output = pd.DataFrame()
    for spin_label, residues in groups:
        column_name = 'avg_dist_{}'.format(spin_label)
        residues = residues.set_index('resid')
        output[column_name] = residues['distance']
        output['restype'] = residues['restype']
    return output


def nmr_to_meld(resid):
    '''
    Convert from NMR to MELD indexing.

    Restype		NMR Index		MELD Index		Note
    ALA			1				0               First protein residue
    ASP         2               1
    GLN         3               2
    ...
    ALA         147             146
    LYS         148             147
    CYS         149             148             Last protein residue

    CA                          149             Calcium
    CA                          150             Calcium
    CA                          151             Calcium
    CA                          152             Calcium

    ALA         201             153             First peptide residue
    ARG         202             154
    ARG         203             155
    ...
    LEU         218             170
    SER         219             171
    SER         220             172             Last peptide residue
    '''
    resid = int(resid)

    if resid <= 149:
        return resid - 1
    elif resid >= 201 and resid <= 220:
        return resid - 48
    else:
        raise RuntimeError('Cannot convert resid {} from NMR to MELD'.format(resid))


def nmr_to_sim_pdb(resid):
    '''
    Convert from NMR to simulated PDB indexing.

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
    missing_residues = [1, 2, 3, 4, 147, 148, 149, 201]
    if resid in missing_residues:
        return None
    elif resid >= 5 and resid <= 146:
        return resid - 4
    elif resid >= 202 and resid <= 220:
        return resid - 55
    else:
        raise RuntimeError('Cannot convert resid {} from NMR to Simulated PDB'.format(resid))


def main():
    # Load the excel file
    excel = pd.ExcelFile('CaCaM2smMLCK_06102015.xls')

    # read in relaxation rates
    df = load_r2(excel)
    df = df.set_index('resid')

    with open('pre_restraints.dat', 'w') as outfile, open('close_residues.tcl', 'w') as tclfile:
        # Loop over the data sets and compare the distances predicted from the NMR data
        # to the distances from the PDB file.
        ds_names = ['S17C', 'T34C', 'N42C', 'N53C', 'R86C', 'T110C', 'T117C', 'E127C', 'Q143C', 'C149']
        for ds_name in ds_names:
            # figure out the name of our new columns
            ond_resid = re.sub(r'[^\d]', '', ds_name)

            # Load in the PRE dataset and merge into the data frame
            pre = load_data(excel, ds_name)
            pre = pre.set_index('resid')
            pre = df.merge(pre, how='outer', left_index=True, right_index=True)

            # compute the average and standard deviations of the distances
            # based on gaussain noise in the peak intensities
            para = pre['paramagnetic']
            dia = pre['diamagnetic']
            r2 = pre['r2']

            dists = compute_distances(para, dia, r2)
            pre['dist'] = dists

            short_pre = pre[pre['dist'] <= 1.5]

            # output the restraints to a file
            for index, row in short_pre.iterrows():
                fixed_index = nmr_to_meld(index)
                fixed_ond_index = nmr_to_meld(ond_resid)
                print('{}\tOND\t{}\tN\t 1.6 250'.format(fixed_ond_index + 1, fixed_index + 1), file=outfile)

            # output vmd commands to add bonds
            # skip 149 because its not in the pdb file
            for index, row in short_pre.iterrows():
                fixed_index = nmr_to_sim_pdb(index)
                fixed_ond_index = nmr_to_sim_pdb(ond_resid)
                if fixed_index is None:
                    print('Could not find {}'.format(index))
                    continue
                if fixed_ond_index is None:
                    print('Could not find {}'.format(ond_resid))
                    continue

                print(
                    'set sel1 [atomselect top "resid {} and name OND"]'.format(fixed_ond_index),
                    file=tclfile)
                print(
                    'set sel2 [atomselect top "resid {} and name N"]'.format(fixed_index),
                    file=tclfile)
                print('label add Bonds 0/[$sel1 list] 0/[$sel2 list]', file=tclfile)


if __name__ == '__main__':
    main()
