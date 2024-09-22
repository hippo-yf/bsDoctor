
from traceback import print_tb
from scipy.stats import gaussian_kde
from src.config import params
from src.utils import *


# DNAme (DP>=k) vs depth
def plot_me_vs_depth_ge_k() -> None:
    dict_genome_me = params['dict_genome_me']
    dict_genome_covnC = params['dict_genome_covnC']
    DP = params['MAXDP_IN_FIG']
    img_dir = params['img_dir']
    save_svg = params['save_svg']

    for cg in CONTEXTS:
        fig, ax = plt.subplots(figsize=(5, 3))
        x = np.arange(1, DP)
        yd = nandivide(dict_genome_me[cg]['double'][:DP], dict_genome_covnC[cg]['double'][:DP])
        yw = nandivide(dict_genome_me[cg]['W'][:DP], dict_genome_covnC[cg]['W'][:DP])
        yc = nandivide(dict_genome_me[cg]['C'][:DP], dict_genome_covnC[cg]['C'][:DP])

        ax.plot(x, yw[1:], '.-', c=COLS[1], markersize=5, label='Watson strand')
        ax.plot(x, yc[1:], '.-', c=COLS[0], markersize=5, label='Crick strand')
        ax.plot(x, yd[1:], '.-', c=COL_gray, markersize=5, label='double strands')
        ax.legend()
        ax.set_xlabel('depth threshold')
        ax.set_ylabel('methylation level')
        filename = f'{img_dir}/meth-{cg}-vs-dp-threshold'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        if save_svg:
            plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    return None

# DNAme (DP=k) vs depth
def plot_me_vs_depth_eq_k() -> None:
    dict_genome_me = params['dict_genome_me']
    dict_genome_covnC = params['dict_genome_covnC']
    DP = params['MAXDP_IN_FIG']
    img_dir = params['img_dir']

    for cg in CONTEXTS:
        fig, ax = plt.subplots(figsize=(5, 3))
        x = np.arange(1, DP)
        yd = nandivide(depthDiff(dict_genome_me[cg]['double'][:DP]), 
                    depthDiff(dict_genome_covnC[cg]['double'][:DP])
                    )
        yw = nandivide(depthDiff(dict_genome_me[cg]['W'][:DP]), 
                    depthDiff(dict_genome_covnC[cg]['W'][:DP])
                    )
        yc = nandivide(depthDiff(dict_genome_me[cg]['C'][:DP]), 
                    depthDiff(dict_genome_covnC[cg]['C'][:DP])
                    )

        ax.plot(x, yw[1:], '.-', c=COLS[1], markersize=5, label='Watson strand')
        ax.plot(x, yc[1:], '.-', c=COLS[0], markersize=5, label='Crick strand')
        ax.plot(x, yd[1:], '.-', c=COL_gray, markersize=5, label='double strands')
        ax.legend()
        ax.set_xlabel('depth')
        ax.set_ylabel('methylation level')
        filename = f'{img_dir}/meth-{cg}-at-dp-k'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        if params['save_svg']:
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

        ax.plot(xw[1:], yw[1:], '.-', c=COLS[1], markersize=5, label='Watson strand')
        ax.plot(xc[1:], yc[1:], '.-', c=COLS[0], markersize=5, label='Crick strand')
        ax.plot(xd[1:], yd[1:], '.-', c=COL_gray, markersize=5, label='double strands')
        ax.legend()
        ax.set_xlabel('missing rate')
        ax.set_ylabel('methylation level')
        filename = f'{img_dir}/meth-{cg}-vs-missing-rate'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    return None

def plot_me_and_covrate_vs_cytosine_density() -> None:
    # dict_binning = params['dict_binning']
    dict_bin_me = params['dict_bin_me']
    dict_bin_covn = params['dict_bin_covn']
    dict_bin_cden = params['dict_bin_cden']
    img_dir = params['img_dir']
    chrs = params['chrs_valid']
    # MAX_DP_BY_FIG = params['MAX_DP_BY_FIG']

    for cg in CONTEXTS:
        for strand in STRANDS:
            me_all = []
            covn_all = []
            cden_all = []
            for chr in chrs:
                me_all.append(dict_bin_me[cg][strand][chr])
                covn_all.append(dict_bin_covn[cg][strand][chr])
                cden_all.append(dict_bin_cden[cg][strand][chr])
            me_all = np.vstack(me_all)
            covn_all = np.vstack(covn_all)
            cden_all = np.hstack(cden_all)
            
            for dp in range(1, 11):
                i = covn_all[:,dp] >= 5
                density = cden_all[i] * 1000
                meth = me_all[i,dp]
                covrate = covn_all[i,dp]/covn_all[i,0]
                density2 = density
                # sampling 3000 points
                if len(meth) > 4000:
                    j = GEN.choice(len(meth), size=3000, replace=False)
                    density = density[j]
                    meth = meth[j]

                ## meth vs density
                xy = np.vstack([density, meth])
                # print(density, meth)
                d = gaussian_kde(xy)(xy)
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

                filename = f'{img_dir}/meth-vs-{cg}-density-of-{strand}-strand-dp-ge{dp}'
                plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
                if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
                plt.close()

                ## cov vs density
                if len(covrate) > 4000:
                    j = GEN.choice(len(covrate), size=3000, replace=False)
                    density2 = density2[j]
                    covrate = covrate[j]
                xy = np.vstack([density2, covrate])
                d = gaussian_kde(xy)(xy)
                idx = d.argsort()
                x, y, d = density2[idx], covrate[idx], d[idx]

                fig, ax = plt.subplots(figsize=(5, 4))
                plt.scatter(x, y, c=d, s=1, cmap='Spectral_r')
                # plt.ylim(0, 1)
                plt.xlim(np.min(density2)-1, np.max(density2)+1)
                plt.colorbar()
                plt.xlabel('cytosine density')
                plt.ylabel('cytosine coverage rate')

                filename = f'{img_dir}/covrate-vs-{cg}-density-of-{strand}-strand-dp-ge{dp}'
                plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
                if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
                plt.close()
    return None
