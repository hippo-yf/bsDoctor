
import scipy
import numpy as np
from src.config import params
from src.utils import *

def plot_depth_vs_cytosine_density():

    dict_binning = params['dict_binning']
    img_dir = params['img_dir']

    gen = np.random.Generator(np.random.PCG64())

    for cg in CONTEXTS:
        for strand in STRANDS:
            density = []
            depth = []
            for key, value in dict_binning.items():
                value = dict_binning[key]
                dp = {
                    'double': value.dp, 'W': value.dpW, 'C': value.dpC
                }
                nC = {
                    'CG': {'double': value.nCG, 'W': value.nCGW, 'C': value.nCGC},
                    'CHG': {'double': value.nCHG, 'W': value.nCHGW, 'C': value.nCHGC},
                    'CHH': {'double': value.nCHH, 'W': value.nCHHW, 'C': value.nCHHC},
                }
                if value.length >= 100:
                    den = nC[cg][strand] / value.length
                    density.append(den)
                    depth.append(dp[strand] / value.length)
            density = np.asarray(density) * 1000 # density per 1000 bp
            depth = np.asarray(depth)

            # sampling 3000 points
            if len(depth) > 3000:
                i = gen.choice(len(depth), size=3000, replace=False)
                density = density[i]
                depth = depth[i]

            ## depth vs density
            xy = np.vstack([density, depth])
            d = scipy.stats.gaussian_kde(xy)(xy)
            idx = d.argsort()
            x, y, d = density[idx], depth[idx], d[idx]

            fig, ax = plt.subplots(figsize=(5, 4))
            plt.scatter(x, y, c=d, s=1, cmap='Spectral_r')
            plt.ylim(0, np.quantile(y, 0.995))
            plt.xlim(np.min(density)-1, np.quantile(x, 0.995))

            plt.colorbar()
            plt.xlabel('cytosine density')
            plt.ylabel('mean depth')
            
            filename = f'{img_dir}/depth-vs-{cg}-density-of-{strand}-strand'
            plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
            plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
            plt.close()

def plot_depth_dist_by_low_high_me():
    dict_binning = params['dict_binning']
    img_dir = params['img_dir']
    DP = params['MAXDP_IN_FIG']

    x = 0
    y = 0
    y2 = 0
    ym = 0 # mean
    nCGLowMeth = 0 
    nCGHighMeth = 0

    for key, value in dict_binning.items():
        # if chr != key[0]: continue
        value = dict_binning[key]
        nCGLowMeth += value.nCGLowMeth
        nCGHighMeth += value.nCGHighMeth

    x = np.arange(DP) + 1
    y = nCGLowMeth/np.sum(nCGLowMeth)
    y2 = nCGHighMeth/np.sum(nCGHighMeth)

    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(x, y[:DP], '.-', c=COLS[0], linewidth=1, markersize=5, label='hypo-meth CGs (<0.3)')
    ax.plot(x, y2[:DP], '.-', c=COLS[1], linewidth=1, markersize=5, label='hyper-meth CGs (>0.7)')
    ax.legend()
    plt.xlabel('methylation read depth')
    plt.ylabel('proportion')

    filename = f'{img_dir}/depth-of-low-high-meth-CGs'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()

def plot_covrate_vs_depth_of_cytosine():
    dict_genome_covnC = params['dict_genome_covnC']
    dict_Cs = params['dict_Cs']
    img_dir = params['img_dir']
    DP_xdepth = params['MAXDP_IN_FIG']

    for cg in CONTEXTS:
        x = np.arange(DP_xdepth) + 1
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.plot(x, dict_genome_covnC[cg]['W'][:DP_xdepth]/dict_Cs[cg]['W'], '.-', c=COLS[1], alpha=1, linewidth=1, markersize=5, label='Watson strand')
        ax.plot(x, dict_genome_covnC[cg]['C'][:DP_xdepth]/dict_Cs[cg]['C'], '.-', c=COLS[0], alpha=1, linewidth=1, markersize=5, label='Crick strand')
        ax.plot(x, dict_genome_covnC[cg]['double'][:DP_xdepth]/dict_Cs[cg]['double'], '.-', c=COL_gray, linewidth=1, markersize=5, label='double strands')
        ax.legend()
        plt.xlabel('methylation read depth')
        plt.ylabel('genome coverage')

        filename = f'{img_dir}/genome-{cg}-coverage-vs-depth'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')

    plt.close()

def plot_depth_dist_of_cytosine():
    dict_genome_covnC = params['dict_genome_covnC']
    img_dir = params['img_dir']
    DP_xdepth = params['MAXDP_IN_FIG']

    for cg in CONTEXTS:
        x = np.arange(DP_xdepth) + 1
        fig, ax = plt.subplots(figsize=(5, 3))
        ax.plot(x, depthDiff(dict_genome_covnC[cg]['W'])[:DP_xdepth] / dict_genome_covnC[cg]['W'][0], '.-', c=COLS[1], alpha=1, linewidth=1, markersize=5, label='Watson strand')
        ax.plot(x, depthDiff(dict_genome_covnC[cg]['C'])[:DP_xdepth] / dict_genome_covnC[cg]['C'][0], '.-', c=COLS[0], alpha=1, linewidth=1, markersize=5, label='Crick strand')
        ax.plot(x, depthDiff(dict_genome_covnC[cg]['double'])[:DP_xdepth] / dict_genome_covnC[cg]['double'][0], '.-', c=COL_gray, linewidth=1, markersize=5, label='double strands')
        ax.legend()
        plt.xlabel('methylation read depth')
        plt.ylabel('proportion')

        filename = f'{img_dir}/genome-{cg}-depth-distribution'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()

def plot_covrate_vs_depth_of_whole_genome():
    genome_covW = params['genome_covW']
    genome_covC = params['genome_covC']
    genome_cov = params['genome_cov']
    length = params['length']
    DP_xdepth = params['MAXDP_IN_FIG']
    img_dir = params['img_dir']

    x = np.arange(DP_xdepth) + 1
    L = np.sum(length)
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(x, genome_covW[:DP_xdepth]/L, '.-', c=COLS[1], alpha=1, linewidth=1, markersize=5, label='Watson strand')
    ax.plot(x, genome_covC[:DP_xdepth]/L, '.-', c=COLS[0], alpha=1, linewidth=1, markersize=5, label='Crick strand')
    ax.plot(x, genome_cov[:DP_xdepth]/L, '.-', c=COL_gray, linewidth=1, markersize=5, label='double strands')
    ax.legend()
    plt.xlabel('depth threshold')
    plt.ylabel('genome coverage')

    filename = f'{img_dir}/whole-genome-coverage-vs-depth'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()

def plot_depth_dist_of_whole_genome():
    genome_covW = params['genome_covW']
    genome_covC = params['genome_covC']
    genome_cov = params['genome_cov']
    DP_xdepth = params['MAXDP_IN_FIG']
    img_dir = params['img_dir']

    x = np.arange(DP_xdepth) + 1
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.plot(x, depthDiff(genome_covW)[:DP_xdepth]/genome_covW[0], '.-', c=COLS[1], alpha=1, linewidth=1, markersize=5, label='Watson strand')
    ax.plot(x, depthDiff(genome_covC)[:DP_xdepth]/genome_covC[0], '.-', c=COLS[0], alpha=1, linewidth=1, markersize=5, label='Crick strand')
    ax.plot(x, depthDiff(genome_cov)[:DP_xdepth]/genome_cov[0], '.-', c=COL_gray, linewidth=1, markersize=5, label='double strands')
    ax.legend()
    plt.ylim(bottom=0)
    plt.xlabel('depth')
    plt.ylabel('proportion')

    filename = f'{img_dir}/whole-genome-depth-distribution'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()

def plot_depth_watson_vs_crick():
    dict_binning = params['dict_binning']
    img_dir = params['img_dir']

    depthW = []
    depthC = []
    for key, value in dict_binning.items():
        value = dict_binning[key]
        if value.length >= 20:
            depthW.append(value.dpW/value.length)
            depthC.append(value.dpC/value.length)
    depthW = np.asarray(depthW)
    depthC = np.asarray(depthC)

    # sampling 3000 points
    if len(depthW) > 3000:
        i = GEN.choice(len(depthW), size=3000, replace=False)
        depthW = depthW[i]
        depthC = depthC[i]

    ## depth vs density
    xy = np.vstack([depthW, depthC])
    d = scipy.stats.gaussian_kde(xy)(xy)
    idx = d.argsort()
    x, y, d = depthW[idx], depthC[idx], d[idx]

    fig, ax = plt.subplots(figsize=(5, 4))
    plt.scatter(x, y, c=d, s=1, cmap='Spectral_r')
    abline(0, 1, color=COL_gray)
    plt.ylim(0, np.quantile(y, 0.99))
    plt.xlim(0, np.quantile(x, 0.99))
    plt.colorbar()
    plt.xlabel('mean depth of Watson strand')
    plt.ylabel('mean depth of Crick strand')

    filename = f'{img_dir}/depth-of-two-strands'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()

def plot_depth_overall_vs_me():
    dict_binning = params['dict_binning']
    img_dir = params['img_dir']

    depthW = []
    depthC = []
    # depthMe = []
    depthMeW = []
    depthMeC = []
    for key, value in dict_binning.items():
        value = dict_binning[key]
        if value.length >= 20:
            depthW.append(value.dpW/value.length)
            depthC.append(value.dpC/value.length)
            # depthMe.append((value.dpCG+value.dpCHG+value.dpCHH)/(value.nCG+value.nCHG+value.nCHH))
            depthMeW.append((value.dpCGW+value.dpCHGW+value.dpCHHW)/(value.nCGW+value.nCHGW+value.nCHHW))
            depthMeC.append((value.dpCGC+value.dpCHGC+value.dpCHHC)/(value.nCGC+value.nCHGW+value.nCHHC))
    depthW = np.asarray(depthW)
    depthC = np.asarray(depthC)
    # depthMe = np.asarray(depthMe)
    depthMeW = np.asarray(depthMeW)
    depthMeC = np.asarray(depthMeC)

    # sampling 3000 points
    if len(depthW) > 3000:
        i = GEN.choice(len(depthW), size=3000, replace=False)
        depthW = depthW[i]
        depthC = depthC[i]
        depthMeW = depthMeW[i]
        depthMeC = depthMeC[i]

    #### Watson
    xy = np.vstack([depthW, depthMeW])
    d = scipy.stats.gaussian_kde(xy)(xy)
    idx = d.argsort()
    x, y, d = depthMeW[idx], depthW[idx], d[idx]

    fig, ax = plt.subplots(figsize=(5, 4))
    plt.scatter(x, y, c=d, s=1, cmap='Spectral_r')
    abline(0, 1, color=COL_gray)
    lim = max(np.quantile(y, 0.995), np.quantile(x, 0.995))
    plt.ylim(0, lim)
    plt.xlim(0, lim)
    plt.colorbar()
    plt.ylabel('mean depth of Watson strand')
    plt.xlabel('mean methylation depth of Watson strand')
    
    filename = f'{img_dir}/depth-vs-meth-depth-of-watson-strand'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()

    #### Crick
    xy = np.vstack([depthC, depthMeC])
    d = scipy.stats.gaussian_kde(xy)(xy)
    idx = d.argsort()
    x, y, d = depthMeC[idx], depthC[idx], d[idx]

    fig, ax = plt.subplots(figsize=(5, 4))
    plt.scatter(x, y, c=d, s=1, cmap='Spectral_r')
    abline(0, 1, color=COL_gray)
    lim = max(np.quantile(y, 0.995), np.quantile(x, 0.995))
    plt.ylim(0, lim)
    plt.xlim(0, lim)
    plt.colorbar()
    plt.ylabel('mean depth of Crick strand')
    plt.xlabel('mean methylation depth of Crick strand')

    filename = f'{img_dir}/depth-vs-meth-depth-of-crick-strand'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
