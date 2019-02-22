#!/usr/bin/env python


import pandas as pd
import numpy as np
import scipy.optimize as opt
import math
import re
import scipy.spatial.distance as sd
import random
from bokeh.plotting import figure, output_file, show
from bokeh.models.sources import ColumnDataSource
from bokeh.models import  HoverTool, Label, Span, Range1d
from bokeh.layouts import gridplot
from matplotlib import pyplot as pp


def load_crystal_distances():
    AVG_FILE = '../ExtractPreDistances/AverageDistances/average_distances.dat'
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


def load_nmr_data(excel):
    df = pd.read_excel(excel, 'CaM2smMLCK', converters={'Res No.': lambda x: int(x)})

    # rename the columns
    df = df.rename(columns={
        'Res No.': 'resid',
        'Name': 'res_name',
        'PRE': 'gamma_S17',
        'flag': 'flag_S17',
        'PRE.1': 'gamma_T34',
        'flag.1': 'flag_T34',
        'PRE.2': 'gamma_N42',
        'flag.2': 'flag_N42',
        'PRE.3': 'gamma_N53',
        'flag.3': 'flag_N53',
        'PRE.4': 'gamma_R86',
        'flag.4': 'flag_R86',
        'PRE.5': 'gamma_T110',
        'flag.5': 'flag_T110',
        'PRE.6': 'gamma_T117',
        'flag.6': 'flag_T117',
        'PRE.7': 'gamma_E127',
        'flag.7': 'flag_E127',
        'PRE.8': 'gamma_Q143',
        'flag.8': 'flag_Q143',
        'PRE.9': 'gamma_C149',
        'flag.9': 'flag_C149'})

    # drop columns we don't care about
    df = df[['resid', 'res_name', 'flag_common',
             'gamma_S17', 'flag_S17',
             'gamma_T34', 'flag_T34',
             'gamma_N42', 'flag_N42',
             'gamma_N53', 'flag_N53',
             'gamma_R86', 'flag_R86',
             'gamma_T110', 'flag_T110',
             'gamma_T117', 'flag_T117',
             'gamma_E127', 'flag_E127',
             'gamma_Q143', 'flag_Q143',
             'gamma_C149', 'flag_C149']]

    # throw out blank rows and residues that have overlaps
    df = df[df['flag_common'] == 1]
    df = df.set_index('resid')

    return df


def get_distances(gamma, flag):
    k = 1.23e-32            # cm**6 s**-2
                            # from Cordina et al, Biochemistry 2013, 52, 1950-1962.
    tau_c = 9.5e-9          # 9.5 ns
    omega = 600e6           # 600 MHz
    gamma = gamma

    r = (k/gamma * (4 * tau_c + 3 * tau_c / (1 + omega**2 * tau_c**2))) ** (1.0 / 6.0)
    r = r / 1e-7            # convert from cm to nm
    r[flag == 0] = np.NAN
    r[flag == -1] = 1.0
    r[r > 3] = 3.0
    return r


def nmr_to_meld(resid):
    if resid <= 149:
        return resid
    elif resid >= 201 and resid <= 220:
        return resid - 47
    else:
        raise ValueError('Cannot convert resid {}.'.format(resid))


def main():
    # Load the excel file
    excel = pd.ExcelFile('CaM_PRE_Iwahara_method.xls')

    # Compute the distances to paramagnetic centers from the pdb file
    # and return them in a dataframe
    xtal = load_crystal_distances()

    nmr = load_nmr_data(excel)

    # Loop over the data sets and compare the distances predicted from the NMR data
    # to the distances from the PDB file.
    ds_names = ['S17', 'T34', 'N42', 'N53', 'R86', 'T110', 'T117', 'E127', 'Q143', 'C149']
    ond_resids = []
    for ds_name in ds_names:
        # figure out the name of our new columns
        ond_resid = re.sub(r'[^\d]', '', ds_name)
        ond_resids.append(ond_resid)
        gamma_name = 'gamma_{}'.format(ds_name)
        flag_name = 'flag_{}'.format(ds_name)
        r_name = 'r_nmr_{}'.format(ond_resid)

        # compute the distances
        nmr[r_name] = get_distances(nmr[gamma_name], nmr[flag_name])

    df = pd.merge(left=nmr, right=xtal, how='inner', left_index=True, right_index=True)

    # now let's write out a list of all of the restraints
    good = 0
    total = 0
    for ond_id in ond_resids:
        for resid, row in df.iterrows():
            r_nmr = row['r_nmr_{}'.format(ond_id)]
            if ond_id == '149':
                r_cryst = np.nan
            else:
                r_cryst = row['avg_dist_{}'.format(ond_id)]

            # skip over cases where the nmr data is bad
            if np.isnan(r_nmr):
                continue

            # compute the allowed distance range based on r_nmr
            if r_nmr < 1.2:
                r_min = 0.
                r_max = 1.7
            elif r_nmr > 1.2 and r_nmr < 2.0:
                r_min = r_nmr - 0.5
                r_max = r_nmr + 0.5
            else:
                r_min = 1.5
                r_max = 999.

            if not np.isnan(r_cryst):
                total += 1
                if r_cryst > r_min and r_cryst < r_max:
                    good += 1

            print '{}\tOND\t{}\tN\t{:8.3f}\t{:8.3f}\t250.'.format(
                ond_id, nmr_to_meld(resid), r_min, r_max)


    # We now have all of the data loaded in one big dataframe,
    # and we're going to use bokeh to plot it. We'll store the
    # output in plot.html
    output_file('plot.html')

    TOOLS = "tap,help,hover"

    source = ColumnDataSource(data=df)
    plots = []
    for resid, mut_name in zip(ond_resids, ds_names):
        # skip C149 because it's not in the crystal structure
        if mut_name == 'C149':
            continue


        p = figure(plot_width=250, plot_height=250,
                tools=TOOLS)

        p.patch([0, 1.2, 1.2, 0], [0, 0, 1.7, 1.7], color='green', alpha=0.1)
        p.patch([1.2, 2.0, 2.0, 1.2], [0.7, 1.5, 2.5, 1.7], color='green', alpha=0.1)
        p.patch([2.0, 5.0, 5.0, 2.0], [1.5, 1.5, 5.0, 5.0], color='green', alpha=0.1)
        # Draw +/- 0.4 angstrom lines.
        # p.line([0, 4.5], [0.4, 4.9], color='grey')
        # p.line([0, 4.5], [-0.4, 4.1], color='grey')

        # Plot the predicted vs actual distance.
        # The plots will be linked because they all share the same
        # datasource.
        p.circle('r_nmr_{}'.format(resid),
                 'avg_dist_{}'.format(resid),
                 source=source,
                 name='distance')

        # Set the tool-tips
        hover = p.select(dict(type=HoverTool))
        hover.tooltips = [
            ('resid', '@resid'),
            ('restype', '@restype'),
            ('pre', '@r_nmr_{}'.format(resid)),
            ('xtal', '@avg_dist_{}'.format(resid)),
            ('I_para', '@para_{}'.format(resid)),
            ('I_dia', '@dia_{}'.format(resid)),
            ('r2', '@r2')
        ]
        hover.names = ['distance']

        # Add a label
        label = Label(x=0.6, y=4.0, text=mut_name, text_color='grey', text_align='center')
        p.add_layout(label)

        p.x_range = Range1d(0, 3.05)
        p.y_range = Range1d(0, 5.00)

        plots.append(p)

    grid = gridplot(plots, ncols=3)
    show(grid)


if __name__ == '__main__':
    main()
