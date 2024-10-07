
import tqdm
from src.utils import *
from src.config import params, data
from src.coverage import CovContig, GenomicIntervalGenerator
from src.updateBinning import update_contig

def compt_plastid() -> None:
    chr_plastid = params['chr_plastid']
    fa = params['fa']
    bam = params['bam']

    if not chr_plastid in fa.references:
        data['plastid_is_covered'] = 0
        return None
    
    dict_plastid = CovContig(contig=chr_plastid)
    intervals = iter(GenomicIntervalGenerator(
        fa, 
        chrs= chr_plastid, 
        start = 0,    
        end = params['MAX_COORDINATE'],
        step=500_000,
        spacing=0
        ))

    intervals_list = list(intervals)
    for i in tqdm.trange(len(intervals_list), desc='Sampling plastid: '):
        detailedIntvl = bam.detailedCoverageContig(intervals_list[i])
        update_contig(dict_plastid, detailedIntvl)
    
    params['dict_plastid'] = dict_plastid

    # whether plastid DNA is sequenced
    if dict_plastid.covnC < 10:
        data['plastid_is_covered'] = 0
    else:
        data['plastid_is_covered'] = 1

        data['plastid_nCorG'] = fi(dict_plastid.nC)
        data['plastid_covnCorG'] = fi(dict_plastid.covnC)
        data['plastid_covnCorG_prop'] = fp(dict_plastid.covnC/dict_plastid.nC)
        data['plastid_length'] = fi(dict_plastid.length.sum())
        data['plastid_covn'] = fi(dict_plastid.cov.sum())
        data['plastid_cov_prop'] = fp(dict_plastid.cov.sum()/dict_plastid.length.sum())
        # data['plastid_size'] = fi(dict_plastid.length.max())

        i = dict_plastid.length > 0
        plastid_median_dp = float(np.median(dict_plastid.dp[i]/dict_plastid.length[i]))
        data['plastid_median_dp'] = ff(plastid_median_dp)
        # plastid_mean_dp = float(np.mean(dict_plastid.dp[i]/dict_plastid.length[i]))
        plastid_mean_dp = float(np.sum(dict_plastid.dp)/np.sum(dict_plastid.length))
        params['plastid_mean_dp'] = plastid_mean_dp
        data['plastid_mean_dp'] = ff2(plastid_mean_dp, 1)
        
        ## bs rate of plastid DNA
        bs_rate_plastid = -1
        data['bsrate_plastid'] = "nan"
        plastid_me = float(dict_plastid.meC/dict_plastid.covnC/10000)
        data['plastid_me'] = ff(plastid_me)
        # if dict_plastid.covnC > 100: # at least 100 Cs covered
        bs_rate_plastid = 1 - plastid_me
        data['bsrate_plastid'] = fp(bs_rate_plastid)
            
        # base error rate by plastid DNA
        i = dict_plastid.dp20 > 0
        data['err_rate_plastid'] = "nan"
        if sum(i) > 0:
            error_rate_plastid = np.mean(dict_plastid.misbase[i]/dict_plastid.dp20[i])
            data['err_rate_plastid'] = fp(error_rate_plastid)
        data['err_rate_plastid']
    return None

def plot_plastid_depth_binning() -> None:
    dict_plastid = params['dict_plastid']
    binsContig = params['binsContig']
    binSizeContig = params['binSizeContig']
    chr_plastid = params['chr_plastid']
    img_dir = params['img_dir']

    x = np.arange(binsContig[chr_plastid])*binSizeContig[chr_plastid]
    x, prefix = prefixBpSize(x)

    yw = dict_plastid.dpW/dict_plastid.length
    yc = dict_plastid.dpC/dict_plastid.length
    yd = dict_plastid.dp/dict_plastid.length
    ylim = np.quantile(np.block([yw,yc]), 0.99)

    figheight = 3*0.5 + 1
    figwidth = len(x)/200 + 1 + 0.2

    fig, axs = plt.subplots(3, 1, sharex=True, sharey=False, figsize=(figwidth, figheight))
    fig.subplots_adjust(hspace=0)
    plt.xlim(-1, max(x)+1)

    axs[0].plot(x, np.repeat(np.median(yw), len(x)), '--', c='gray')
    axs[0].scatter(x, np.fmin(ylim, yw), c=COLS[1], s=1)
    axs[0].set_ylim(0, ylim)
    axs[0].text(axs[0].get_xlim()[1], axs[0].get_ylim()[1], f'med: {np.median(yw):.0f}', horizontalalignment='right', verticalalignment='top')

    axs[1].plot(x, np.repeat(np.median(yc), len(x)), '--', c='gray')
    axs[1].scatter(x, np.fmin(ylim, yc), c=COLS[0], s=1)
    axs[1].set_ylim(0, ylim)
    axs[1].text(axs[1].get_xlim()[1], axs[1].get_ylim()[1], f'med: {np.median(yc):.0f}', horizontalalignment='right', verticalalignment='top')
    axs[1].set_ylabel('read depth')

    axs[2].plot(x, np.repeat(np.median(yd), len(x)), '--', c='gray')
    axs[2].scatter(x, np.fmin(ylim, yd), c=COL_gray, s=1)
    axs[2].text(axs[2].get_xlim()[1], axs[2].get_ylim()[1], f'med: {np.median(yd):.0f}', horizontalalignment='right', verticalalignment='top')

    plt.xlabel(f'coordinate ({prefix})')

    filename = f'{img_dir}/depth-bin-plastid'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

def plot_plastid_base_error_rate() -> None:
    dict_plastid = params['dict_plastid']
    img_dir = params['img_dir']
    
    fig, ax = plt.subplots(figsize=(5,3))
    if np.sum(dict_plastid.dp20 > 0) > 0:
        i = dict_plastid.dp20 > 0
        ax.hist(dict_plastid.misbase[i]/dict_plastid.dp20[i], bins=21, density=True, color=COLS[0])
        plt.xlabel('bin-wise base error rate')
        plt.ylabel('density')

        filename = f'{img_dir}/base-error-rate-dist-plastid'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    return None
