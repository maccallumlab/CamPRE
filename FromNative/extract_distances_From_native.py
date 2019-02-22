#!/usr/bin/env python


import pandas as pd
import numpy as np
import scipy.optimize as opt
from matplotlib import pyplot as pp
import seaborn as sns
import math


def load_r2(excel):
    df = pd.read_excel(excel, 'S17C', header=None)
    df.columns = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'residue', 'r2', 'j']
    r2 = df.drop(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'j'], 1)
    r2 = r2.dropna()
    r2['residue'] = [int(r) for r in r2['residue']]
    r2 = r2.set_index('residue')
    return r2


def load_data(excel, sheetname):
    df = pd.read_excel(excel, sheetname, header=None)
    df.columns = [
        'residue',
        'diamagnetic',
        'c',
        'paramagnetic',
        'e',
        'f',
        'g',
        'h',
        'i',
        'j']
    # drop 'g' first, because it's all NaN
    ir = df.drop('g', 1)
    # now drop anything left with NaN
    ir = ir.dropna()
    # now drop the columns we don't need
    ir = ir.drop(['c', 'e', 'f', 'h', 'i', 'j'], 1)
    ir['residue'] = [int(r) for r in ir['residue']]
    ir = ir.set_index('residue')
    return ir


def compute_ratio(ir):
    ratios = ir['paramagnetic'].values / ir['diamagnetic'].values
    # normalize so the top 10% of ration have median of 1
    n = len(ratios) // 3
    correction = np.median(sorted(ratios)[-n:])
    pp.plot(ratios)
    pp.axhline(correction)
    pp.show()
    print correction
    ratios = ratios / correction
    ir['ratio'] = ratios
    ir['ratio'][ir['ratio'] > 1.2] = np.nan
    return ir.dropna()


def compute_gamma(ir):
    gammas = []
    for ratio, r2 in zip(ir['ratio'], ir['r2']):
        def func(gamma):
            return (r2 * math.exp(-gamma * 10e-3) / (r2 + gamma) - ratio)**2
        x = opt.newton(func, 4., maxiter=500)
        gammas.append(x)
    ir['gamma'] = gammas
    ir[ir['gamma'] < 0] = 0
    return ir


def compute_distance(pre):
    K = 1.23e-32 * 1e-12 # cm^6 s^-2 * m^6/cm^6 = m^6 s^-2
    tc = 8e-9 # ns
    omega_H = 700e6 # s^-1
    f = 4 * tc + 3 * tc / (1 + (omega_H * tc)**2 )
    rs = []
    g = K / (1.0e-9)**6 * f
    for gamma in pre['gamma']:
        r6 = K * f / gamma
        r = r6**(1.0 / 6.0) * 1e9
        rs.append(r)
    pre['r'] = rs
    return pre


def process_bounds(pre, cutoff1, cutoff2):
    def get_lower(r):
        if r < cutoff1:
            return 0.
        elif r > cutoff2:
            return cutoff2
        else:
            return max(0, r - 0.4)

    def get_upper(r):
        if r < cutoff1:
            return cutoff1
        elif r > cutoff2:
            return 999
        else:
            return r + 0.4

    lower = [get_lower(r) for r in pre['r'].values]
    upper = [get_upper(r) for r in pre['r'].values]
    pre['lower'] = lower
    pre['upper'] = upper
    return pre


excel = pd.ExcelFile('../RawData/CaCaM2smMLCK_06102015.xls')
r2 = load_r2(excel)

DS_NAME = 'N42C'
INDEX = 42
pre = load_data(excel, DS_NAME)
pre = pd.merge(pre, r2, left_index=True, right_index=True)

pre = compute_ratio(pre)
pre = compute_gamma(pre)
pre = compute_distance(pre)
pre = process_bounds(pre, cutoff1=1.2, cutoff2=2.0)

for index, row in pre.iterrows():
    fmt = '{:d}\tOND\t{:d}\tN\t{:.3f}\t{:.3f}'
    print fmt.format(17, index, row['lower'], row['upper'])
