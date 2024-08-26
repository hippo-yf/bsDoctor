
import math
# from os import path
import os
import sys
import gzip
import io
import random

import numpy as np
import matplotlib.pyplot as plt
import pysam
from numpy.typing import NDArray
from typing import NamedTuple, List, Tuple, Dict, Any

# rng = np.random.default_rng()
GEN = np.random.Generator(np.random.PCG64())

COLS = ['#1B98E0', '#df5353', '#00A896', '#ff9c33', '#32a896']
COLS_AREA = ['#49afe9', '#e36969', 'ffa84d']
COL_gray = '#666666'
CONTEXTS = ['CG', 'CHG', 'CHH']
STRANDS = ['double', 'W', 'C']
BASES = ['A', 'T', 'C', 'G']


def axlimit(x):
    k = int(np.log10(x))
    if x < 1: k -= 1
    return (int(x * 10**(1-k))+1) * 10**(k-1)

def prefixBpSize(x):
    M = np.max(x)
    if M >= 1e9: 
        x = x/1e9
        prefix = 'Gbp'
    elif M >= 1e6:
        x = x/1e6
        prefix = 'Mbp'
    elif M >= 1e3:
        x = x/1e3
        prefix = 'Kbp'
    else: prefix = 'bp'
    return x, prefix

# divide keeping nan's
def nandivide(x, y):
    i = (y != 0)
    q = np.zeros(np.shape(x), dtype=float)
    q[i] = x[i]/y[i]
    q[~i] = math.nan
    return q

# trim values larger than quantile 0.99
def trimQuantile(x, q=0.99):
    return np.fmin(x, np.quantile(x, q))

def depthDiff(x: NDArray):
    return -np.diff(np.hstack([x, 0]))
def depthCulSum(x: NDArray):
    return np.cumsum(x[::-1])[::-1]

def signif(x, digit=3):
    return np.floor(np.array(x)*10**digit)/10**digit


# format numbers in strings
def fi(x: int):
    return f'{x:,}'

def ff(x: float, k:int=4):
    s = f'%.{k}g'
    return s % x

def ff2(x: float, k:int=2):
    s = f'%.{k}f'
    return s % x

def fp(x, k:int=4):
    return '{1:.{0:d}g}'.format(k, x*100) # '%.{k}g%%'

def fp2(x, k:int=2):
    return '{1:.{0:d}f}'.format(k, x*100)

def reverseComp(str: str) -> str:
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join([comp[s] for s in str[::-1]])

# 3-nucleotide context, CG/CHG/CHH etc.

CG_CONTEXT_FORWARD_HASH = {'CGA':'CG', 'CGT':'CG', 'CGC':'CG', 'CGG':'CG', 'CGN':'CG', # 4 CG
                           'CAG':'CHG', 'CTG':'CHG', 'CCG':'CHG',          # 3 CHG
                           'CAA':'CHH', 'CAT':'CHH', 'CAC':'CHH',          # 9 CHH
                           'CTA':'CHH', 'CTT':'CHH', 'CTC':'CHH',
                           'CCA':'CHH', 'CCT':'CHH', 'CCC':'CHH'
                           }

CG_CONTEXT_REVERSE_HASH = {'ACG':'CG', 'TCG':'CG', 'CCG':'CG', 'GCG':'CG', 'NCG':'CG', # 4 CG
                           'CAG':'CHG', 'CTG':'CHG', 'CGG':'CHG',          # 3 CHG
                           'AAG':'CHH', 'ATG':'CHH', 'AGG':'CHH',          # 9 CHH
                           'TAG':'CHH', 'TTG':'CHH', 'TGG':'CHH',
                           'GAG':'CHH', 'GTG':'CHH', 'GGG':'CHH'
                           }

# 2-nucleotide context of reverse strand

DI_CONTEXT_REVERSE_HASH = {'AG':'CT', 'TG':'CA', 'CG':'CG', 'GG':'CC'}

def as_bool(x: str):
    x = x.upper()
    if x in ["TRUE", "T", "YES", "Y"] : return True
    if x in ["FALSE", "F", "NO", "N"] : return False
    return None

def myOpenFile(file: str):
    if file == '':
        return None
    if file == '-':
        outfile = sys.stdout
    elif file.endswith('.gz'):
        outfile = gzip.open(file, 'wt')
    else:
        outfile = io.open(file, 'wt')
    return outfile

def myMakeDirs(dir: str):
    if not os.path.exists(dir):
        os.makedirs(dir)


def getBins(length, bins):
    size = math.ceil(length/bins)
    if size > 10: size = math.ceil(size/10)*10
    bins = math.ceil(length/size)
    return (size, bins)

def sampleGenes(l: list, k:int = 1000):
    k = min(k, len(l))
    s = random.sample(range(len(l)), k)
    s.sort()
    return [l[i] for i in s]

# DNAme adjustment for incomplete bs conversion
def adjust_me(me, bsrate):
    return 1 - (1 - me)/bsrate

def relative_freq_CpG_motif (x, df=4):
    # df: number of free bases
    return x/np.sum(x) * 4**df

## utils for figures
def abline(intercept, slope, **kw):
    """Plot a line from slope and intercept"""
    axes = plt.gca()
    x_vals = np.array(axes.get_xlim())
    y_vals = intercept + slope * x_vals
    plt.plot(x_vals, y_vals, '--', **kw)


