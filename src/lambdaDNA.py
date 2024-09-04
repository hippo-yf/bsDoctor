
import tqdm
from src.utils import *
from src.config import params, data
from src.coverage import CovLambda, GenomicIntervalGenerator
from src.updateBinning import update_lambda

def compt_lambda() -> None:
    chr_lambda = params['chr_lambda']
    fa = params['fa']
    bam = params['bam']

    dict_lambda = CovLambda(contig=chr_lambda)
    intervals = iter(GenomicIntervalGenerator(
        fa, 
        chrs= chr_lambda, 
        start = 0,    
        end = params['MAX_COORDINATE'],
        step=5_000,
        spacing=0
        ))

    intervals_list = list(intervals)
    for i in tqdm.trange(len(intervals_list), desc='Sampling control DNA: '):
        detailedIntvl = bam.detailedCoverageContig(intervals_list[i])
        update_lambda(dict_lambda, detailedIntvl)
    
    params['dict_lambda'] = dict_lambda

    # whether lambda DNA is sequenced
    if dict_lambda.cov.sum() < 10:
        data['lambda_is_covered'] = 0
    else:
        data['lambda_is_covered'] = 1

        data['lambda_nCorG'] = fi(dict_lambda.nC)
        data['lambda_covnCorG'] = fi(dict_lambda.covnC)
        data['lambda_covnCorG_prop'] = fp(dict_lambda.covnC/dict_lambda.nC)
        data['lambda_length'] = fi(dict_lambda.length.sum())
        data['lambda_covn'] = fi(dict_lambda.cov.sum())
        data['lambda_cov_prop'] = fp(dict_lambda.cov.sum()/dict_lambda.length.sum())
        data['lambda_size'] = fi(dict_lambda.length.max())

        i = dict_lambda.length > 0
        lambda_median_dp = np.median(dict_lambda.dp[i]/dict_lambda.length[i])
        data['lambda_median_dp'] = ff(float(lambda_median_dp))
        # lambda_mean_dp = np.mean(dict_lambda.dp[i]/dict_lambda.length[i])
        lambda_mean_dp = np.sum(dict_lambda.dp)/np.sum(dict_lambda.length)
        data['lambda_mean_dp'] = ff(float(lambda_mean_dp))
        
        ## bs rate of lambda DNA
        bs_rate_lambda = -1
        data['bsrate_lambda'] = "nan"
        if dict_lambda.covnC > 100: # at least 100 Cs covered
            bs_rate_lambda = 1 - dict_lambda.meC/dict_lambda.covnC
            params['bs_rate_lambda'] = bs_rate_lambda # for last choice of conversion rate
            data['bsrate_lambda'] = fp(bs_rate_lambda)
            
        # base error rate by lambda DNA
        i = dict_lambda.dp20 > 0
        data['err_rate_lambda'] = "nan"
        if sum(i) > 0:
            error_rate_lambda = np.mean(dict_lambda.misbase[i]/dict_lambda.dp20[i])
            data['err_rate_lambda'] = fp(error_rate_lambda)
        data['err_rate_lambda']
    return None

def plot_lambda_depth_binning() -> None:
    dict_lambda = params['dict_lambda']
    binsContig = params['binsContig']
    binSizeContig = params['binSizeContig']
    chr_lambda = params['chr_lambda']
    img_dir = params['img_dir']
    save_svg = params['save_svg']

    x = np.arange(binsContig[chr_lambda])*binSizeContig[chr_lambda]
    x, prefix = prefixBpSize(x)

    yw = dict_lambda.dpW/dict_lambda.length
    yc = dict_lambda.dpC/dict_lambda.length
    yd = dict_lambda.dp/dict_lambda.length
    ylim = np.quantile(np.block([yw,yc]), 0.99)

    figheight = 3*0.5 + 1
    figwidth = len(x)/200 + 1 + 0.2

    fig, axs = plt.subplots(3, 1, sharex=True, sharey=False, figsize=(figwidth, figheight))
    fig.subplots_adjust(hspace=0)
    plt.xlim(-1, max(x)+1)

    axs[0].plot(x, np.repeat(np.median(yw), len(x)), '--', c='gray')
    axs[0].scatter(x, np.fmin(ylim, yw), c=COLS[1], s=1)
    axs[0].set_ylim(0, ylim)
    axs[0].text(2, axs[0].get_ylim()[1], f'med: {np.median(yw):.0f}', horizontalalignment='left', verticalalignment='top')

    axs[1].plot(x, np.repeat(np.median(yc), len(x)), '--', c='gray')
    axs[1].scatter(x, np.fmin(ylim, yc), c=COLS[0], s=1)
    axs[1].set_ylim(0, ylim)
    axs[1].text(2, axs[1].get_ylim()[1], f'med: {np.median(yc):.0f}', horizontalalignment='left', verticalalignment='top')
    axs[1].set_ylabel('read depth')

    axs[2].plot(x, np.repeat(np.median(yd), len(x)), '--', c='gray')
    axs[2].scatter(x, np.fmin(ylim, yd), c=COL_gray, s=1)
    axs[2].text(2, axs[2].get_ylim()[1], f'med: {np.median(yd):.0f}', horizontalalignment='left', verticalalignment='top')

    plt.xlabel(f'coordinate ({prefix})')

    filename = f'{img_dir}/depth-bin-lambda'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    if save_svg:
        plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

def plot_lambda_base_error_rate() -> None:
    dict_lambda = params['dict_lambda']
    img_dir = params['img_dir']
    save_svg = params['save_svg']
    
    fig, ax = plt.subplots(figsize=(5,3))
    if np.sum(dict_lambda.dp20 > 0) > 0:
        ax.hist(nandivide(dict_lambda.misbase, dict_lambda.dp20), bins=21, density=True, color=COLS[0])
        plt.xlabel('bin-wise base error rate')
        plt.ylabel('density')

        filename = f'{img_dir}/base-error-rate-dist-lambda'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        if save_svg:
            plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    return None
    
