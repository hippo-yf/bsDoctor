
import tqdm
from src.utils import *
from src.config import params, data
from src.coverage import CovContig, GenomicIntervalGenerator
from src.updateBinning import update_contig

def compt_MT() -> None:
    chr_MT = params['chr_MT']
    fa = params['fa']
    bam = params['bam']

    dict_MT = CovContig(contig=chr_MT)
    intervals = iter(GenomicIntervalGenerator(
        fa, 
        chrs=chr_MT, 
        start=0,
        end = params['MAX_COORDINATE'],
        step=5_000,
        spacing=0
        ))

    intervals_list = list(intervals)
    for i in tqdm.trange(len(intervals_list), desc='Sampling MT: '):
        detailedIntvl = bam.detailedCoverageContig(intervals_list[i])
        update_contig(dict_MT, detailedIntvl)

    params['dict_MT'] = dict_MT

    # whether lambda DNA is sequenced
    if dict_MT.cov.sum() < 10:
        data['mt_is_covered'] = 0
    else:
        data['mt_is_covered'] = 1

    data['mt_nCorG'] = fi(dict_MT.nC)
    data['mt_covnCorG'] = fi(dict_MT.covnC)
    data['mt_covnCorG_prop'] = fp(dict_MT.covnC/dict_MT.nC)
    data['mt_length'] = fi(dict_MT.length.sum())
    data['mt_covn'] = fi(dict_MT.cov.sum())
    data['mt_cov_prop'] = fp(dict_MT.cov.sum()/dict_MT.length.sum())
    data['mt_bin_size'] = fi(dict_MT.length.max())
    
    i = dict_MT.length > 0
    mt_median_dp = np.median(dict_MT.dp[i]/dict_MT.length[i])
    data['mt_median_dp'] = ff(float(mt_median_dp))
    # mt_mean_dp = np.mean(dict_MT.dp[i]/dict_MT.length[i])
    mt_mean_dp = np.sum(dict_MT.dp)/np.sum(dict_MT.length)
    params['mt_mean_dp'] = mt_mean_dp
    data['mt_mean_dp'] = ff2(float(mt_mean_dp), 1)

    ## conversion rate by MT
    mt_me = float(dict_MT.meC/dict_MT.covnC/10000)
    data['mt_me'] = ff(mt_me)
    bs_rate_MT = 1 - mt_me
    params['bs_rate_mt'] = bs_rate_MT
    # dp_MT = dict_MT.dp/dict_MT.length
    data['bsrate_mt'] = fp(bs_rate_MT)

    ## eror rate by MT
    i = dict_MT.dp20 > 0
    data['err_rate_mt'] = 'nan'
    if sum(i) > 0:
        error_rate_MT = np.mean(dict_MT.misbase[i]/dict_MT.dp20[i])
        data['err_rate_mt'] = fp(error_rate_MT)
    data['err_rate_mt']
    return None

def plot_mt_depth_binning() -> None:
    dict_MT = params['dict_MT']
    binsContig = params['binsContig']
    binSizeContig = params['binSizeContig']
    chr_MT = params['chr_MT']
    img_dir = params['img_dir']

    x = np.arange(binsContig[chr_MT])*binSizeContig[chr_MT]
    x, prefix = prefixBpSize(x)

    yw = dict_MT.dpW/dict_MT.length
    yc = dict_MT.dpC/dict_MT.length
    yd = dict_MT.dp/dict_MT.length
    ylim = np.quantile(np.block([yw,yc]), 0.99)

    figheight = 3*0.5 + 1
    figwidth = len(x)/200 + 1 + 0.2

    fig, axs = plt.subplots(3, 1, sharex=True, sharey=False, figsize=(figwidth, figheight))
    fig.subplots_adjust(hspace=0)

    axs[0].plot(x, np.repeat(np.median(yw), len(x)), '--', c='gray')
    axs[0].scatter(x, np.fmin(ylim, yw), c=COLS[1], s=1)
    axs[0].set_xlim(0, max(x))
    axs[0].set_ylim(0, ylim)
    axs[0].text(axs[0].get_xlim()[1], axs[0].get_ylim()[1], f'med: {np.median(yw):.0f}', horizontalalignment='right', verticalalignment='top')

    axs[1].plot(x, np.repeat(np.median(yc), len(x)), '--', c='gray')
    axs[1].scatter(x, np.fmin(ylim, yc), c=COLS[0], s=1)
    axs[1].set_xlim(0, max(x))
    axs[1].set_ylim(0, ylim)
    axs[1].text(axs[1].get_xlim()[1], axs[1].get_ylim()[1], f'med: {np.median(yc):.0f}', horizontalalignment='right', verticalalignment='top')
    axs[1].set_ylabel(f'read depth')

    axs[2].plot(x, np.repeat(np.median(yd), len(x)), '--', c='gray')
    axs[2].scatter(x, np.fmin(ylim, yd), c=COL_gray, s=1)
    axs[2].set_xlim(0, max(x))
    axs[2].text(axs[2].get_xlim()[1], axs[2].get_ylim()[1], f'med: {np.median(yd):.0f}', horizontalalignment='right', verticalalignment='top')

    plt.xlabel(f'coordinate ({prefix})')

    filename = f'{img_dir}/depth-bin-MT'
    # filename = f'img/depth-bin-MT'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

def plot_mt_base_error_rate() -> None:
    dict_MT = params['dict_MT']
    img_dir = params['img_dir']

    fig, ax = plt.subplots(figsize=(5,3))
    if np.sum(dict_MT.dp20 > 0) > 10:
        i = dict_MT.dp20 > 0
        ax.hist(dict_MT.misbase[i]/dict_MT.dp20[i], bins=21, density=True, color=COLS[0])
        plt.xlabel('bin-wise base error rate')
        plt.ylabel('density')

        filename = f'{img_dir}/base-error-rate-dist-MT'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    return None
