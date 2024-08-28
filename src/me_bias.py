
import scipy
from src.config import params
from src.utils import *


# DNAme (DP>=k) vs depth
def plot_me_vs_depth_eq_k() -> None:
    dict_genome_me = params['dict_genome_me']
    dict_genome_covnC = params['dict_genome_covnC']
    DP = params['MAXDP_IN_FIG']
    img_dir = params['img_dir']

    for cg in CONTEXTS:
        fig, ax = plt.subplots(figsize=(5, 3))
        x = np.arange(DP) + 1
        yd = nandivide(dict_genome_me[cg]['double'][:DP], dict_genome_covnC[cg]['double'][:DP])
        yw = nandivide(dict_genome_me[cg]['W'][:DP], dict_genome_covnC[cg]['W'][:DP])
        yc = nandivide(dict_genome_me[cg]['C'][:DP], dict_genome_covnC[cg]['C'][:DP])

        ax.plot(x, yw, '.-', c=COLS[1], markersize=5, label='Watson strand')
        ax.plot(x, yc, '.-', c=COLS[0], markersize=5, label='Crick strand')
        ax.plot(x, yd, '.-', c=COL_gray, markersize=5, label='double strands')
        ax.legend()
        ax.set_xlabel('depth threshold')
        ax.set_ylabel('methylation level')
        filename = f'{img_dir}/meth-{cg}-vs-dp-threshold'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    return None

# DNAme (DP=k) vs depth
def plot_me_vs_depth_ge_k() -> None:
    dict_genome_me = params['dict_genome_me']
    dict_genome_covnC = params['dict_genome_covnC']
    DP = params['MAXDP_IN_FIG']
    img_dir = params['img_dir']

    for cg in CONTEXTS:
        fig, ax = plt.subplots(figsize=(5, 3))
        x = np.arange(DP) + 1
        yd = nandivide(depthDiff(dict_genome_me[cg]['double'][:DP]), 
                    depthDiff(dict_genome_covnC[cg]['double'][:DP])
                    )
        yw = nandivide(depthDiff(dict_genome_me[cg]['W'][:DP]), 
                    depthDiff(dict_genome_covnC[cg]['W'][:DP])
                    )
        yc = nandivide(depthDiff(dict_genome_me[cg]['C'][:DP]), 
                    depthDiff(dict_genome_covnC[cg]['C'][:DP])
                    )

        ax.plot(x, yw, '.-', c=COLS[1], markersize=5, label='Watson strand')
        ax.plot(x, yc, '.-', c=COLS[0], markersize=5, label='Crick strand')
        ax.plot(x, yd, '.-', c=COL_gray, markersize=5, label='double strands')
        ax.legend()
        ax.set_xlabel('depth')
        ax.set_ylabel('methylation level')
        filename = f'{img_dir}/meth-{cg}-at-dp-k'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    return None

def plot_me_vs_missing() -> None:
    dict_genome_me = params['dict_genome_me']
    dict_genome_covnC = params['dict_genome_covnC']
    dict_Cs = params['dict_Cs']
    DP = params['MAXDP_IN_FIG']
    img_dir = params['img_dir']

    for cg in CONTEXTS:
        fig, ax = plt.subplots(figsize=(5, 3))
        # dp = np.arange(DP) + 1
        xd = 1 - dict_genome_covnC[cg]['double'][:DP] / dict_Cs[cg]['double']
        yd = nandivide(dict_genome_me[cg]['double'][:DP], dict_genome_covnC[cg]['double'][:DP])
        xw = 1 - dict_genome_covnC[cg]['W'][:DP] / dict_Cs[cg]['W']
        yw = nandivide(dict_genome_me[cg]['W'][:DP], dict_genome_covnC[cg]['W'][:DP])
        xc = 1 - dict_genome_covnC[cg]['C'][:DP] / dict_Cs[cg]['C']
        yc = nandivide(dict_genome_me[cg]['C'][:DP], dict_genome_covnC[cg]['C'][:DP])

        ax.plot(xw, yw, '.-', c=COLS[1], markersize=5, label='Watson strand')
        ax.plot(xc, yc, '.-', c=COLS[0], markersize=5, label='Crick strand')
        ax.plot(xd, yd, '.-', c=COL_gray, markersize=5, label='double strands')
        ax.legend()
        ax.set_xlabel('missing rate')
        ax.set_ylabel('methylation level')
        filename = f'{img_dir}/meth-{cg}-vs-missing-rate'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    return None

def plot_me_and_covrate_vs_cytosine_density() -> None:
    dict_binning = params['dict_binning']
    img_dir = params['img_dir']

    gen = np.random.Generator(np.random.PCG64())
    for cg in CONTEXTS:
        for strand in STRANDS:
            for dp in range(10):
                density = []
                meth = []
                density2 = [] # for cov 
                covrate = []
                for key, value in dict_binning.items():
                    value = dict_binning[key]
                    me = {
                        'CG': {'double': value.meCG[dp], 'W': value.meCGW[dp], 'C': value.meCGC[dp]},
                        'CHG': {'double': value.meCHG[dp], 'W': value.meCHGW[dp], 'C': value.meCHGC[dp]},
                        'CHH': {'double': value.meCHH[dp], 'W': value.meCHHW[dp], 'C': value.meCHHC[dp]}
                    }
                    nC = {
                        'CG': {'double': value.nCG, 'W': value.nCGW, 'C': value.nCGC},
                        'CHG': {'double': value.nCHG, 'W': value.nCHGW, 'C': value.nCHGC},
                        'CHH': {'double': value.nCHH, 'W': value.nCHHW, 'C': value.nCHHC},
                    }
                    covnC = {
                        'CG': {'double': value.covnCG[dp], 'W': value.covnCGW[dp], 'C': value.covnCGC[dp]},
                        'CHG': {'double': value.covnCHG[dp], 'W': value.covnCHGW[dp], 'C': value.covnCHGC[dp]},
                        'CHH': {'double': value.covnCHH[dp], 'W': value.covnCHHW[dp], 'C': value.covnCHHC[dp]},
                    }
                    if value.length >= 100 and value.covnCGW[dp] >= 10 and value.covnCGC[dp] >= 10 and value.covnCHGW[dp] >= 20 and value.covnCHGC[dp] >= 20 and value.covnCHHW[dp] >= 30 and value.covnCHHC[dp] >= 30:
                        density.append(nC[cg][strand] / value.length)
                        meth.append(me[cg][strand] / covnC[cg][strand])
                    if value.length >= 200 and value.nCG >= 10 and value.nCHG >= 20 and value.nCHH >= 30:
                        density2.append(nC[cg][strand] / value.length)
                        covrate.append(covnC[cg][strand] / nC[cg][strand])

                density = np.asarray(density) * 1000 # density per 1000 bp
                meth = np.asarray(meth) # meth
                density2 = np.asarray(density2) * 1000
                covrate = np.asarray(covrate) # cov rate
                # sampling 2000 points
                if len(meth) > 3000:
                    i = gen.choice(len(meth), size=3000, replace=False)
                    density = density[i]
                    meth = meth[i]

                ## meth vs density
                xy = np.vstack([density, meth])
                d = scipy.stats.gaussian_kde(xy)(xy)
                # Sort the points by density
                # so that the densest points are plotted at last
                idx = d.argsort()
                x, y, d = density[idx], meth[idx], d[idx]

                fig, ax = plt.subplots(figsize=(5, 4))
                plt.scatter(x, y, c=d, s=1, cmap='Spectral_r')
                # plt.ylim(0, 1)
                plt.xlim(np.min(density)-1, np.max(density)+1)
                plt.colorbar()
                plt.xlabel('cytosine density')
                plt.ylabel('methylation level')

                filename = f'{img_dir}/meth-vs-{cg}-density-of-{strand}-strand-dp-ge{dp+1}'
                plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
                plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
                plt.close()

                ## cov vs density
                if len(covrate) > 3000:
                    i = gen.choice(len(covrate), size=3000, replace=False)
                    density2 = density2[i]
                    covrate = covrate[i]
                xy = np.vstack([density2, covrate])
                d = scipy.stats.gaussian_kde(xy)(xy)
                idx = d.argsort()
                x, y, d = density2[idx], covrate[idx], d[idx]

                fig, ax = plt.subplots(figsize=(5, 4))
                plt.scatter(x, y, c=d, s=1, cmap='Spectral_r')
                # plt.ylim(0, 1)
                plt.xlim(np.min(density2)-1, np.max(density2)+1)
                plt.colorbar()
                plt.xlabel('cytosine density')
                plt.ylabel('cytosine coverage rate')

                filename = f'{img_dir}/covrate-vs-{cg}-density-of-{strand}-strand-dp-ge{dp+1}'
                plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
                plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
                plt.close()
    return None
