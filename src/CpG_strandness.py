
from matplotlib.pyplot import plot
from src.config import params, data
from src.utils import *



def plot_heatmap_stranded_CpG_depth():
    stranded_CG_depth = params['stranded_CG_depth']
    img_dir = params['img_dir']

    # n = stranded_CG_depth.shape[0]
    x = np.arange(5) *20

    fig, ax = plt.subplots()
    stranded_CG_depth[0,0] = 0
    im = ax.imshow(np.log10(1+stranded_CG_depth), cmap='Spectral_r')
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    ax.set_xticks(x, labels=x)
    ax.set_yticks(x, labels=x)
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('log10 (#CpGs + 1)', rotation=-90, va="bottom")
    plt.xlabel('methylation read depth of Watson strand')
    plt.ylabel('methylation read depth of Crick strand')

    filename = f'{img_dir}/heatmap-cg-stranded-depth'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()

def plot_bar_CpG_depth_difference():
    stranded_CG_depth = params['stranded_CG_depth']
    img_dir = params['img_dir']
    DP_xdepth = params['MAXDP_IN_FIG']

    cov_diff = np.zeros((199,), dtype=int)
    shape = np.shape(stranded_CG_depth)
    for i in range(shape[0]):
        for j in range(shape[1]):
            if (i,j) == (0,0): continue
            cov_diff[i-j+99] += stranded_CG_depth[i,j]
    dp_diff = np.arange(-99, 100)
    cov_diff = cov_diff/cov_diff.sum()

    dp = min(DP_xdepth+20, 100)

    fig, ax = plt.subplots(figsize=(7, 2.5))
    ax.bar(dp_diff[(100-dp):(100+dp)], cov_diff[(100-dp):(100+dp)], color=COLS[0])
    plt.xlabel('depth difference of a CpG')
    plt.ylabel('proportion')

    filename = f'{img_dir}/cg-stranded-depth-difference'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()

    #### stats of |depth diff|
    # mean/sd, M(mean)AD/sd of cov diff
    a = np.inner(dp_diff, cov_diff)
    sda = math.sqrt(np.inner(dp_diff**2, cov_diff) - a**2)
    b = np.inner(np.abs(dp_diff), cov_diff)
    sdb = math.sqrt(np.inner(dp_diff**2, cov_diff) - b**2)
    params['mean_cpg_depth_diff'] = a
    params['sd_cpg_depth_diff'] = sda
    params['mean_abscpg_depth_diff'] = b
    params['sd_abs_cpg_depth_diff'] = sdb

def plot_bar_double_srtanded_cpg():
    stranded_CG_depth = params['stranded_CG_depth']
    img_dir = params['img_dir']
    DP_xdepth = params['MAXDP_IN_FIG']

    prop_double_cov = np.zeros((DP_xdepth,))
    shape = np.shape(stranded_CG_depth)
    for i in np.arange(1, DP_xdepth+1):
        b = stranded_CG_depth[i:, i:].sum() # both
        a = stranded_CG_depth[i:,:].sum() + stranded_CG_depth[:,i:].sum() - b
        prop_double_cov[i-1] = b/a

    fig, ax = plt.subplots(figsize=(5,3))
    ax.bar(np.arange(1, DP_xdepth+1), prop_double_cov, color=COLS[0])
    plt.xlabel('min{Watson depth, Crick depth}')
    plt.ylabel('proportion of double-stranded CpGs')

    filename = f'{img_dir}/cg-prop-of-double-stranded-coverage'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()

def plot_heatmap_stranded_CpG_meth():
    stranded_CG_meth = params['stranded_CG_meth']
    img_dir = params['img_dir']

    # DP >= 1
    # n = stranded_CG_meth.shape[0]
    # n = params['meth_bins']
    n = 20
    x = np.arange(n)
    labels = [
        f'[{k/n:.2f}, {(k+1)/n:.2f})' 
        if k<n-1 else f'[{k/n:.2f}, {(k+1)/n:.2f}]' 
        for k in np.arange(n)
        ]
    fig, ax = plt.subplots()
    im = ax.imshow(np.log10(1+stranded_CG_meth), cmap='Spectral_r')
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    ax.set_xticks(x, labels=labels)
    ax.set_yticks(x, labels=labels)
    plt.xlabel('DNAme of Crick strand')
    plt.ylabel('DNAme of Watson strand')

    ax.tick_params(axis='x', labelrotation= 90)
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('log10 (#CpGs + 1)', rotation=-90, va="bottom")

    filename = f'{img_dir}/heatmap-cg-stranded-meth'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()

def plot_heatmap_stranded_meth_diff():
    CGmeth_diff_by_depth = params['CGmeth_diff_by_depth']
    img_dir = params['img_dir']
    DP_xdepth = params['MAXDP_IN_FIG']

    # n = stranded_CG_meth.shape[0]
    n = 20
    # n = params['meth_bins']

    y = np.arange(2, 40, step=3)
    dp = min(100, DP_xdepth+20)
    x = np.arange(start=0, stop=dp-1, step=5)
    labels = [
        f'[{k/n:5.2f}, {(k+1)/n:5.2f})' 
        if k<n-1 else f'[{k/n:5.2f}, {(k+1)/n:5.2f}]' 
        for k in np.arange(-18, 20, step=3)
        ]
    fig, ax = plt.subplots()
    im = ax.imshow(np.log10(1+CGmeth_diff_by_depth[:,1:dp]), cmap='Spectral_r')
    # im = ax.imshow(np.log10(1+stranded_CG_meth), cmap='Spectral_r')
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    ax.set_yticks(y, labels=labels)
    ax.set_xticks(x, labels=x+1)
    plt.xlabel('min{Watson depth, Crick depth}')
    plt.ylabel('DNAme difference between strands')

    # ax.tick_params(axis='x', labelrotation= 90)
    cbar = ax.figure.colorbar(im, ax=ax, orientation='horizontal', fraction=0.05, pad=0.08)
    # cbar.ax.set_ylabel('log10 (#CpGs + 1)', rotation=-90, va="bottom")
    cbar.ax.set_xlabel('log10 (#CpGs + 1)')

    filename = f'{img_dir}/heatmap-cg-stranded-meth-diff'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()

