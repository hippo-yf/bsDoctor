
from numpy import float64, int64
from src.utils import *
from src.config import params, data

from src.utils import adjust_me

def compt_whole_genome() -> None:
    dict_binning = params['dict_binning']
    MAXDEPTH = params['MAXDEPTH']

    length = []
    totalDepth = []
    totalDepthW = []
    totalDepthC = []
    binning_depth = []
    binning_depthW = []
    binning_depthC = []
    binning_cov = []
    binning_covW = []
    binning_covC = []

    stranded_CG_depth = 0
    stranded_CG_meth = 0
    CGmeth_diff_by_depth = 0

    # whole-genome counts
    meCG: NDArray = np.zeros((MAXDEPTH,), dtype=float64)
    meCGW: NDArray = np.zeros((MAXDEPTH,), dtype=float64)
    meCGC: NDArray = np.zeros((MAXDEPTH,), dtype=float64)
    meCHG: NDArray = np.zeros((MAXDEPTH,), dtype=float64)
    meCHGW: NDArray = np.zeros((MAXDEPTH,), dtype=float64)
    meCHGC: NDArray = np.zeros((MAXDEPTH,), dtype=float64)
    meCHH: NDArray = np.zeros((MAXDEPTH,), dtype=float64)
    meCHHW: NDArray = np.zeros((MAXDEPTH,), dtype=float64)
    meCHHC: NDArray = np.zeros((MAXDEPTH,), dtype=float64)
    covnCG: NDArray = np.zeros((MAXDEPTH,), dtype=int64)
    covnCGW: NDArray = np.zeros((MAXDEPTH,), dtype=int64)
    covnCGC: NDArray = np.zeros((MAXDEPTH,), dtype=int64)
    covnCHG: NDArray = np.zeros((MAXDEPTH,), dtype=int64)
    covnCHGW: NDArray = np.zeros((MAXDEPTH,), dtype=int64)
    covnCHGC: NDArray = np.zeros((MAXDEPTH,), dtype=int64)
    covnCHH: NDArray = np.zeros((MAXDEPTH,), dtype=int64)
    covnCHHW: NDArray = np.zeros((MAXDEPTH,), dtype=int64)
    covnCHHC: NDArray = np.zeros((MAXDEPTH,), dtype=int64)

    # in total
    nCG = 0
    nCGW = 0
    nCGC = 0
    nCHG = 0
    nCHGW = 0
    nCHGC = 0
    nCHH = 0
    nCHHW = 0
    nCHHC = 0
    ATdp = []
    misbase = []
    nCGMethBin = 0
    nCHGMethBin = 0
    nCHHMethBin = 0

    for key, value in dict_binning.items():
        length.append(value.length)
        totalDepth.append(value.dp)
        totalDepthW.append(value.dpW)
        totalDepthC.append(value.dpC)
        binning_depth.append(-np.diff(np.block([value.cov, 0])))
        binning_depthW.append(-np.diff(np.block([value.covW, 0])))
        binning_depthC.append(-np.diff(np.block([value.covC, 0])))
        binning_cov.append(value.cov)
        binning_covW.append(value.covW)
        binning_covC.append(value.covC)
        stranded_CG_depth += value.stranded_CG_depth
        stranded_CG_meth += value.stranded_CG_meth
        CGmeth_diff_by_depth += value.CGmeth_diff_by_depth
        
        nCG += value.nCG
        nCGW += value.nCGW
        nCGC += value.nCGC
        nCHG += value.nCHG
        nCHGW += value.nCHGW
        nCHGC += value.nCHGC
        nCHH += value.nCHH
        nCHHW += value.nCHHW
        nCHHC += value.nCHHC
        covnCG += value.covnCG
        covnCGW += value.covnCGW
        covnCGC += value.covnCGC
        covnCHG += value.covnCHG
        covnCHGW += value.covnCHGW
        covnCHGC += value.covnCHGC
        covnCHH += value.covnCHH
        covnCHHW += value.covnCHHW
        covnCHHC += value.covnCHHC
        meCG += value.meCG
        meCGW += value.meCGW
        meCGC += value.meCGC
        meCHG += value.meCHG
        meCHGW += value.meCHGW
        meCHGC += value.meCHGC
        meCHH += value.meCHH
        meCHHW += value.meCHHW
        meCHHC += value.meCHHC
        ATdp.append(value.ATdp)
        misbase.append(value.misbase)
        nCGMethBin += value.nCGMethBin
        nCHGMethBin += value.nCHGMethBin
        nCHHMethBin += value.nCHHMethBin

    length = np.asarray(length)
    params['length'] = length # length of bins
    binning_depth = np.asarray(binning_depth)
    binning_depthW = np.asarray(binning_depthW)
    binning_depthC = np.asarray(binning_depthC)
    binning_cov = np.asarray(binning_cov)
    binning_covW = np.asarray(binning_covW)
    binning_covC = np.asarray(binning_covC)
    params['binning_covW'] = binning_covW
    params['binning_covC'] = binning_covC
    params['binning_depthW'] = binning_depthW
    params['binning_depthC'] = binning_depthC

    # whole-genome cov
    genome_cov = np.sum(binning_cov, axis=0)
    genome_covW = np.sum(binning_covW, axis=0)
    genome_covC = np.sum(binning_covC, axis=0)
    params['genome_cov'] = genome_cov
    params['genome_covW'] = genome_covW
    params['genome_covC'] = genome_covC

    totalDepth = np.asarray(totalDepth)
    totalDepthW = np.asarray(totalDepthW)
    totalDepthC = np.asarray(totalDepthC)
    params['totalDepth'] = totalDepth
    params['totalDepthW'] = totalDepthW
    params['totalDepthC'] = totalDepthC
    totalBases = totalDepth.sum()
    totalBasesW = totalDepthW.sum()
    totalBasesC = totalDepthC.sum()
    params['totalBases'] = totalBases
    params['totalBasesW'] = totalBasesW
    params['totalBasesC'] = totalBasesC
    ATdp =  np.asarray(ATdp)
    misbase = np.asarray(misbase)
    
    # define DP threshold
    prop_CG = covnCG/covnCG[0]
    DP_xdepth = min(MAXDEPTH, max(10, MAXDEPTH -1 - np.argmax(prop_CG[::-1] >= 0.03)))
    params['MAXDP_IN_FIG'] = DP_xdepth

    #### summary in dict
    # culmulative covered Cs
    dict_genome_covnC = { 
    'CG': {'double': covnCG, 'W': covnCGW, 'C': covnCGC},
    'CHG': {'double': covnCHG, 'W': covnCHGW, 'C': covnCHGC},
    'CHH': {'double': covnCHH, 'W': covnCHHW, 'C': covnCHHC}
    }

    # culmulative me
    dict_genome_me = { 
    'CG': {'double': meCG, 'W': meCGW, 'C': meCGC},
    'CHG': {'double': meCHG, 'W': meCHGW, 'C': meCHGC},
    'CHH': {'double': meCHH, 'W': meCHHW, 'C': meCHHC}
    }

    # genome Cs
    dict_Cs = {
        'CG':{'double': nCG, 'W': nCGW, 'C': nCGC},
        'CHG':{'double': nCHG, 'W': nCHGW, 'C': nCHGC},
        'CHH':{'double': nCHH, 'W': nCHHW, 'C': nCHHC}
    }

    params['dict_genome_covnC'] = dict_genome_covnC
    params['dict_genome_me'] = dict_genome_me
    params['dict_Cs'] = dict_Cs

    # counts of Cs by meth level binning
    methDist = {'CG': nCGMethBin, 'CHG': nCHGMethBin, 'CHH': nCHHMethBin}
    params['methDist'] = methDist

    # CpG strandness
    params['stranded_CG_depth'] = stranded_CG_depth
    params['stranded_CG_meth'] = stranded_CG_meth
    params['CGmeth_diff_by_depth'] = CGmeth_diff_by_depth

    ################################################
    #### for report
    ################################################

    data['nCHH'] = fi(nCHH)
    #### bs rate by CHH weighted by dp
    wdp = np.sqrt(np.arange(MAXDEPTH-1))
    w = -np.diff(covnCHH)*wdp
    mechh = np.sum(-np.diff(meCHH)*wdp)/w.sum()
    bs_rate_chh = 1 - mechh
    data['bsrate_chh'] = fp(bs_rate_chh)
    params['bs_rate_chh'] = bs_rate_chh

    #### adopted bs conversion rate
    # if 'bs_rate_lambda' in params and 'bs_rate_MT' in params:

    if bool(data['include_lambda']) and bool(data['lambda_is_covered']):
        bs_rate = params['bs_rate_lambda']
    elif bool(data['include_mt']) and bool(data['MT_is_covered']):
        bs_rate = params['bs_rate_MT']
    else:
        bs_rate = bs_rate_chh

    params['bs_rate'] = bs_rate
    data['bs_rate'] = fp(bs_rate)

    #### adjust whole-genome DNAme
    data['me_CG'] = fp(meCG[0]/covnCG[0])
    data['me_CG_adj'] = fp(adjust_me(meCG[0]/covnCG[0], bs_rate))
    data['me_CHG'] = fp(meCHG[0]/covnCHG[0])
    data['me_CHG_adj'] = fp(adjust_me(meCHG[0]/covnCHG[0], bs_rate))
    data['me_CHH'] = fp(meCHH[0]/covnCHH[0])
    data['me_CHH_adj'] = fp(adjust_me(meCHH[0]/covnCHH[0], bs_rate))
    # data['bsrate_warning'] = int(max(bs_rate_chh, bs_rate_MT, bs_rate_lambda) < 0.9)
    data['bsrate_warning'] = int(bs_rate < 0.9)
    
    #### depth for all bases
    i = length > 0
    median_dp = np.median(totalDepth[i]/length[i])
    data['median_dp'] = ff(float(median_dp))
    # mean_dp = np.mean(totalDepth[i]/length[i])
    mean_dp = np.sum(totalDepth)/np.sum(length)
    data['mean_dp'] = ff(float(mean_dp))

    #### MT DNA copy number
    if 'mt_mean_dp' in params:
        data['mt_copynum'] = ff(params['mt_mean_dp'] / mean_dp *2, 4)

    #### lambda DNA copy number
    if 'lambda_mean_dp' in params:
        data['lambda_copynum'] = ff(params['lambda_mean_dp'] / mean_dp *2, 4)

    #### plastmid DNA copy number
    if 'plastmid_mean_dp' in params:
        data['plastmid_copynum'] = ff(params['plastmid_mean_dp'] / mean_dp *2, 4)

    #### error rate by A/T sites
    # bins with total AT depth >= 100 
    i = ATdp >= 100
    err_AT = misbase[i]/ATdp[i]
    err_rate_AT = np.mean(err_AT)
    data['err_rate_AT'] = fp(err_rate_AT)
    params['err_AT'] = err_AT
    return None

def plot_base_error_rate_by_AT() -> None:
    err_AT = params['err_AT']
    img_dir = params['img_dir']

    fig, ax = plt.subplots(figsize=(5,3))
    if len(err_AT) > 0:
        ax.hist(err_AT, bins=21, density=True, color=COLS[0])
        plt.xlabel('base error rate')
        plt.ylabel('density')
        filename = f'{img_dir}/base-error-rate-by-AT'
        # filename = f'img/base-error-rate-by-AT'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    return None

def plot_theroretical_me_bias() -> None:
    # DNAme bias
    bs_rate = params['bs_rate']

    rate = np.linspace(0.8, 1, 200)
    meth = np.linspace(0, 1, 200)
    meth0 = 1 - np.outer(1-meth, 1/rate)
    meth_bias = np.outer(meth, np.ones(200)) - meth0

    # constrant: meth + bsrate >= 1
    tmp = np.outer(meth, np.ones(200)) + np.outer(np.ones(200), rate)
    meth_bias[tmp<1] = np.nan

    fig, ax = plt.subplots()
    im = ax.imshow(meth_bias, cmap='Spectral_r')

    # vertical line of bsrate
    if 0.8 < bs_rate < 1:
        # y0 = np.linspace(1, 1-bsrate, num=10)
        y0 = [int(max(5, (1-bs_rate)*200)), 195]
        x0 = [int(bs_rate*200), int(bs_rate*200)]
        plt.plot(x0, y0, color='gray', linestyle='dashed')

    ax.tick_params(top=False, bottom=True, labeltop=False, labelbottom=True)
    x = np.linspace(0, 200-1, 6)
    lx = [f'{x:.2f}' for x in np.linspace(0.8, 1, 6)]
    y = np.linspace(0, 200-1, 6)
    ly = [f'{y:.2f}' for y in np.linspace(0, 1, 6)]
    ax.set_xticks(x, labels=lx)
    ax.set_yticks(y, labels=ly)
    plt.xlabel('Conversion rate')
    plt.ylabel('Unadjusted DNAme level')

    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('DNAme level bias', rotation=-90, va="bottom")
    
    filename = f'{params['img_dir']}/DNAme-bias'
    # filename = f'img/DNAme-bias'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

# count/desity of Cs by meth level binning
def plot_hist_me() -> None:
    methDist = params['methDist']
    img_dir = params['img_dir']

    for key, value in methDist.items():
        count = np.cumsum(value[:,::-1], axis=1)[:,::-1]
        prop = count / np.sum(count, axis=0)
        # ylim = axlimit(np.max(prop))
        shape = np.shape(prop)
        x = np.arange(shape[0])
        x2 = np.arange(0,shape[0]+1, 4)
        xlabels=x2/shape[0]
        for DP in range(shape[1]):
        # for DP in range(1):
            fig, ax = plt.subplots(figsize=(5,2))
            ax.bar(x, prop[:,DP], width=1, align='edge', color='#1B98E0')
            ax.set_xticks(x2, xlabels)
            # ax.set_ylim(0, ylim)
            ax.set_xlabel('mean DNAme level')
            ax.set_ylabel('proportion')

            filename = f'{img_dir}/meth-dist-genome-{key}-dp{DP+1}'
            plt.savefig(filename+'.png', transparent=True, dpi=120, bbox_inches='tight')
            plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
            plt.close()
    return None
