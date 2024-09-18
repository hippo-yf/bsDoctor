
from src.utils import *
from src.config import params, data


def compt_whole_genome() -> None:
    dict_binning = params['dict_binning']
    ploidy = params['ploidy']
    DP = params['MAXDEPTH']
    MAX_DP_BY_FIG = params['MAX_DP_BY_FIG'] + 1
    L = len(dict_binning)

    # length = []
    length: NDArray[int32] = np.zeros((L,), dtype=int32)
    totalDepth: NDArray[int32] = np.zeros((L,), dtype=int32)
    totalDepthW: NDArray[int32] = np.zeros((L,), dtype=int32)
    totalDepthC: NDArray[int32] = np.zeros((L,), dtype=int32)
    # binning_depth = np.zeros((L, DP), dtype=int32)
    # binning_depthW = np.zeros((L, DP), dtype=int32)
    # binning_depthC = np.zeros((L, DP), dtype=int32)
    binning_cov: NDArray[int32] = np.zeros((L, DP), dtype=int32)
    binning_covW: NDArray[int32] = np.zeros((L, DP), dtype=int32)
    binning_covC: NDArray[int32] = np.zeros((L, DP), dtype=int32)

    stranded_CG_depth: NDArray[int64] = np.zeros((100,100), dtype=int64)
    stranded_CG_meth: NDArray[int64] = np.zeros((20,20), dtype=int64)
    CGmeth_diff_by_depth: NDArray[int64] = np.zeros((40,100), dtype=int64)

    # whole-genome counts
    # integral me (X10000)
    imeCG: NDArray[int64] = np.zeros((DP,), dtype=int64)
    imeCGW: NDArray[int64] = np.zeros((DP,), dtype=int64)
    imeCGC: NDArray[int64] = np.zeros((DP,), dtype=int64)
    imeCHG: NDArray[int64] = np.zeros((DP,), dtype=int64)
    imeCHGW: NDArray[int64] = np.zeros((DP,), dtype=int64)
    imeCHGC: NDArray[int64] = np.zeros((DP,), dtype=int64)
    imeCHH: NDArray[int64] = np.zeros((DP,), dtype=int64)
    imeCHHW: NDArray[int64] = np.zeros((DP,), dtype=int64)
    imeCHHC: NDArray[int64] = np.zeros((DP,), dtype=int64)
    covnCG: NDArray[int64] = np.zeros((DP,), dtype=int64)
    covnCGW: NDArray[int64] = np.zeros((DP,), dtype=int64)
    covnCGC: NDArray[int64] = np.zeros((DP,), dtype=int64)
    covnCHG: NDArray[int64] = np.zeros((DP,), dtype=int64)
    covnCHGW: NDArray[int64] = np.zeros((DP,), dtype=int64)
    covnCHGC: NDArray[int64] = np.zeros((DP,), dtype=int64)
    covnCHH: NDArray[int64] = np.zeros((DP,), dtype=int64)
    covnCHHW: NDArray[int64] = np.zeros((DP,), dtype=int64)
    covnCHHC: NDArray[int64] = np.zeros((DP,), dtype=int64)

    # in total
    nCG: int = 0
    nCGW: int = 0
    nCGC: int = 0
    nCHG: int = 0
    nCHGW: int = 0
    nCHGC: int = 0
    nCHH: int = 0
    nCHHW: int = 0
    nCHHC: int = 0
    ATdp: NDArray[int32] = np.zeros((L,), dtype=int32)
    misbase: NDArray[int32] = np.zeros((L,), dtype=int32)
    nCGMethBin: NDArray[int32] = np.zeros((20, MAX_DP_BY_FIG), dtype=int32)
    nCHGMethBin: NDArray[int32] = np.zeros((20, MAX_DP_BY_FIG), dtype=int32)
    nCHHMethBin: NDArray[int32] = np.zeros((20, MAX_DP_BY_FIG), dtype=int32)

    for i, (key, value) in enumerate(dict_binning.items()):
        length[i] = value.length
        totalDepth[i] = value.dpW + value.dpC
        totalDepthW[i] = value.dpW
        totalDepthC[i] = value.dpC

        # binning_depth[i,:] = value.cov
        # binning_depthW[i,:] = value.covW
        # binning_depthC[i,:] = value.covC
        # binning_depth.append(-np.diff(np.block([value.cov, 0])))
        # binning_depthW.append(-np.diff(np.block([value.covW, 0])))
        # binning_depthC.append(-np.diff(np.block([value.covC, 0])))

        #sites covered by dp
        binning_cov[i,:] = value.cov
        binning_covW[i,:] = value.covW
        binning_covC[i,:] = value.covC
        stranded_CG_depth += value.stranded_CG_depth
        stranded_CG_meth += value.stranded_CG_meth
        CGmeth_diff_by_depth += value.CGmeth_diff_by_depth
        
        nCG += value.nCGW + value.nCGC
        nCGW += value.nCGW
        nCGC += value.nCGC
        nCHG += value.nCHGW + value.nCHGC
        nCHGW += value.nCHGW
        nCHGC += value.nCHGC
        nCHH += value.nCHHW + value.nCHHC
        nCHHW += value.nCHHW
        nCHHC += value.nCHHC
        covnCG += value.covnCGW + value.covnCGC
        covnCGW += value.covnCGW
        covnCGC += value.covnCGC
        covnCHG += value.covnCHGW + value.covnCHGC
        covnCHGW += value.covnCHGW
        covnCHGC += value.covnCHGC
        covnCHH += value.covnCHHW + value.covnCHHC
        covnCHHW += value.covnCHHW
        covnCHHC += value.covnCHHC
        imeCG += value.meCGW + value.meCGC
        imeCGW += value.meCGW
        imeCGC += value.meCGC
        imeCHG += value.meCHGW + value.meCHGC
        imeCHGW += value.meCHGW
        imeCHGC += value.meCHGC
        imeCHH += value.meCHHW + value.meCHHC
        imeCHHW += value.meCHHW
        imeCHHC += value.meCHHC
        ATdp[i] = value.ATdp
        misbase[i] = value.misbase
        nCGMethBin += value.nCGMethBin
        nCHGMethBin += value.nCHGMethBin
        nCHHMethBin += value.nCHHMethBin

    # save cummualtive ones 
    covnCG = cumsumrev(covnCG)
    covnCGW = cumsumrev(covnCGW)
    covnCGC = cumsumrev(covnCGC)
    covnCHG = cumsumrev(covnCHG)
    covnCHGW = cumsumrev(covnCHGW)
    covnCHGC = cumsumrev(covnCHGC)
    covnCHH = cumsumrev(covnCHH)
    covnCHHW = cumsumrev(covnCHHW)
    covnCHHC = cumsumrev(covnCHHC)
    # float me
    meCG = cumsumrev(imeCG)/10000
    meCGW = cumsumrev(imeCGW)/10000
    meCGC = cumsumrev(imeCGC)/10000
    meCHG = cumsumrev(imeCHG)/10000
    meCHGW = cumsumrev(imeCHGW)/10000
    meCHGC = cumsumrev(imeCHGC)/10000
    meCHH  = cumsumrev(imeCHH)/10000
    meCHHW = cumsumrev(imeCHHW)/10000
    meCHHC = cumsumrev(imeCHHC)/10000
    # length = np.asarray(length)
    params['length'] = length # length of bins
    # binning_depth = np.asarray(binning_depth)
    # binning_depthW = np.asarray(binning_depthW)
    # binning_depthC = np.asarray(binning_depthC)
    # binning_cov = np.asarray(binning_cov)
    # binning_covW = np.asarray(binning_covW)
    # binning_covC = np.asarray(binning_covC)
    params['binning_cov'] = binning_cov
    params['binning_covW'] = binning_covW
    params['binning_covC'] = binning_covC
    # params['binning_depthW'] = binning_depthW
    # params['binning_depthC'] = binning_depthC

    # whole-genome cov
    genome_cov = np.sum(binning_cov, axis=0)
    genome_covW = np.sum(binning_covW, axis=0)
    genome_covC = np.sum(binning_covC, axis=0)
    params['genome_cov'] = genome_cov
    params['genome_covW'] = genome_covW
    params['genome_covC'] = genome_covC

    params['totalDepth'] = totalDepth
    params['totalDepthW'] = totalDepthW
    params['totalDepthC'] = totalDepthC
    totalBases = totalDepth.sum()
    totalBasesW = totalDepthW.sum()
    totalBasesC = totalDepthC.sum()
    params['totalBases'] = totalBases
    params['totalBasesW'] = totalBasesW
    params['totalBasesC'] = totalBasesC
    # ATdp =  np.asarray(ATdp)
    # misbase = np.asarray(misbase)
    
    ## whole-genome covrate
    # DPs = [0,2,4,9]
    DPs = [1,3,5,10]
    covrate_wg_dp = genome_cov/length.sum()
    data['covrate_wg_dp'] = [fp(covrate_wg_dp[k]) for k in DPs]
    covrate_wg_sp = (genome_covW + genome_covC)/length.sum()/2
    data['covrate_wg_sp'] = [fp(covrate_wg_sp[k]) for k in DPs]

    # define DP threshold in figs
    # prop_CG = covnCG/covnCG[0]
    prop_CG = covnCG[1:]/covnCG[1] # prop over dp>=1
    DP_xdepth = min(DP, max(20, DP -1 - np.argmax(prop_CG[::-1] >= 0.03)))
    # max dp in fig plotted by dp
    params['MAXDP_IN_FIG'] = DP_xdepth

    ##  CG covrate
    covrate_CG = covnCG/nCG
    data['covrate_cg'] = [fp(covrate_CG[k]) for k in DPs]

    ##  CHG covrate
    covrate_CHG = covnCHG/nCHG
    data['covrate_chg'] = [fp(covrate_CHG[k]) for k in DPs]

    ##  CHH covrate
    covrate_CHH = covnCHH/nCHH
    data['covrate_chh'] = [fp(covrate_CHH[k]) for k in DPs]

    # CpGs of double-stranded cov
    prop_double_cov = np.zeros((len(DPs),), dtype=float64)
    for k, dp in enumerate(DPs):
        # dp += 1
        b = stranded_CG_depth[dp:, dp:].sum() # both
        a = stranded_CG_depth[dp:,:].sum() + stranded_CG_depth[:,dp:].sum() - b
        prop_double_cov[k] = b/a
    data['covrate_cg_ds'] = [fp(x) for x in prop_double_cov]

    params['stranded_CG_depth'] = stranded_CG_depth

    ## DNAme level
    data['me_cg_ds'] = [fp(x) for x in nandivide(meCG[DPs], covnCG[DPs])]
    data['me_chg_ds'] = [fp(x) for x in nandivide(meCHG[DPs], covnCHG[DPs])]
    data['me_chh_ds'] = [fp(x) for x in nandivide(meCHH[DPs], covnCHH[DPs])]

    #### summary in dict
    # cumulative covered Cs
    dict_genome_covnC = { 
    'CG': {'double': covnCG, 'W': covnCGW, 'C': covnCGC},
    'CHG': {'double': covnCHG, 'W': covnCHGW, 'C': covnCHGC},
    'CHH': {'double': covnCHH, 'W': covnCHHW, 'C': covnCHHC}
    }

    # cumulative me
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
    wdp = np.sqrt(np.arange(DP-1))
    w = -np.diff(covnCHH)*wdp
    mechh = np.sum(-np.diff(meCHH)*wdp)/w.sum()
    bs_rate_chh = 1 - mechh
    data['bsrate_chh'] = fp(bs_rate_chh)
    params['bs_rate_chh'] = bs_rate_chh

    #### adopted bs conversion rate

    if bool(data['include_lambda']) and bool(data['lambda_is_covered']):
        bs_rate = params['bs_rate_lambda']
    elif bool(data['include_mt']) and bool(data['mt_is_covered']):
        bs_rate = params['bs_rate_mt']
    else:
        bs_rate = bs_rate_chh

    params['bs_rate'] = bs_rate
    data['bs_rate'] = fp(bs_rate)

    #### adjust whole-genome DNAme
    data['me_CG'] = fp(meCG[1]/covnCG[1])
    data['me_CG_adj'] = fp(adjust_me(meCG[1]/covnCG[1], bs_rate))
    data['me_CHG'] = fp(meCHG[1]/covnCHG[1])
    data['me_CHG_adj'] = fp(adjust_me(meCHG[1]/covnCHG[1], bs_rate))
    data['me_CHH'] = fp(meCHH[1]/covnCHH[1])
    data['me_CHH_adj'] = fp(adjust_me(meCHH[1]/covnCHH[1], bs_rate))
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
        data['mt_copynum'] = ff(params['mt_mean_dp'] / mean_dp *ploidy, 4)

    #### lambda DNA copy number
    if 'lambda_mean_dp' in params:
        data['lambda_copynum'] = ff(params['lambda_mean_dp'] / mean_dp *ploidy, 4)

    #### plastmid DNA copy number
    if 'plastmid_mean_dp' in params:
        data['plastmid_copynum'] = ff(params['plastmid_mean_dp'] / mean_dp *ploidy, 4)

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
        if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
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
    if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
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
            if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
            plt.close()
    return None
