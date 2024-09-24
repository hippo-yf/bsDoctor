
import tqdm
from src.utils import *
from src.config import params, data
from src.coverage import GenomicIntervalGeneratorWithinBin, BinCov
from src.updateBinning import update_cgkmer, update_binning_nuclear 


def nuclear_sampling() -> None:
    fa = params['fa']
    bam = params['bam']
    chrs_valid = params['chrs_valid']
    step = params['nuclear_sampling_step']
    binSize = params['binSize']
    spacing = params['nuclear_sampling_spacing']
    include_motif = data['include_motif']

    # intervals = iter(GenomicIntervalGenerator(
    intervals = iter(GenomicIntervalGeneratorWithinBin(
        fa, 
        chrs=chrs_valid,
        start=0,    
        end=params['MAX_COORDINATE'],
        step=step, # 1_000
        spacing=spacing, #10_000
        binSize=binSize
        ))

    intervals_list = list(intervals)
    for i in tqdm.trange(len(intervals_list), desc='Sampling nuclear chr: '):
        detailedIntvl = bam.detailedCoverage(intervals_list[i])
        update_binning_nuclear(detailedIntvl)
        if include_motif: 
            update_cgkmer(detailedIntvl)
    return None

def compt_chr_and_bin_wise() -> None:
    # fa = params['fa']
    reference_length = params['reference_length']
    binSize = params['binSize']
    dict_binning = params['dict_binning']
    chrs_valid = params['chrs_valid']
    DP = params['MAX_DP_BY_FIG'] + 1

    ###########################################################
    #### init bin-wise data
    ###########################################################

    dict_bin_meCG = dict()
    dict_bin_meCGW = dict()
    dict_bin_meCGC = dict()
    dict_bin_meCHG = dict()
    dict_bin_meCHGW = dict()
    dict_bin_meCHGC = dict()
    dict_bin_meCHH = dict()
    dict_bin_meCHHW = dict()
    dict_bin_meCHHC = dict()
    dict_bin_covnCG = dict()
    dict_bin_covnCGW = dict()
    dict_bin_covnCGC = dict()
    dict_bin_covnCHG = dict()
    dict_bin_covnCHGW = dict()
    dict_bin_covnCHGC = dict()
    dict_bin_covnCHH = dict()
    dict_bin_covnCHHW = dict()
    dict_bin_covnCHHC = dict()
    dict_bin_dp = dict()
    dict_bin_dpW = dict()
    dict_bin_dpC = dict()
    dict_bin_CGden = dict() # C density
    dict_bin_CGWden = dict()
    dict_bin_CGCden = dict()
    dict_bin_CHGden = dict()
    dict_bin_CHGWden = dict()
    dict_bin_CHGCden = dict()
    dict_bin_CHHden = dict()
    dict_bin_CHHWden = dict()
    dict_bin_CHHCden = dict()
    # DP_valid = 30
    
    ###########################################################
    #### init chr-wise data
    ###########################################################
    dict_chr_meCG = dict()
    dict_chr_meCGW = dict()
    dict_chr_meCGC = dict()
    dict_chr_meCHG = dict()
    dict_chr_meCHGW = dict()
    dict_chr_meCHGC = dict()
    dict_chr_meCHH = dict()
    dict_chr_meCHHW = dict()
    dict_chr_meCHHC = dict()
    dict_chr_mmeCG = dict()
    dict_chr_mmeCGW = dict()
    dict_chr_mmeCGC = dict()
    dict_chr_mmeCHG = dict()
    dict_chr_mmeCHGW = dict()
    dict_chr_mmeCHGC = dict()
    dict_chr_mmeCHH = dict()
    dict_chr_mmeCHHW = dict()
    dict_chr_mmeCHHC = dict()
    dict_chr_covnCG = dict()
    dict_chr_covnCGW = dict()
    dict_chr_covnCGC = dict()
    dict_chr_covnCHG = dict()
    dict_chr_covnCHGW = dict()
    dict_chr_covnCHGC = dict()
    dict_chr_covnCHH = dict()
    dict_chr_covnCHHW = dict()
    dict_chr_covnCHHC = dict()

    # init
    for i, chr in enumerate(chrs_valid):
        nbin = int(np.ceil(reference_length(chr)/binSize))
        shape = (nbin, DP)
        dict_bin_meCG[chr] = np.zeros(shape, dtype=float32)
        dict_bin_meCGW[chr] = np.zeros(shape, dtype=float32)
        dict_bin_meCGC[chr] = np.zeros(shape, dtype=float32)
        dict_bin_meCHG[chr] = np.zeros(shape, dtype=float32)
        dict_bin_meCHGW[chr] = np.zeros(shape, dtype=float32)
        dict_bin_meCHGC[chr] = np.zeros(shape, dtype=float32)
        dict_bin_meCHH[chr] = np.zeros(shape, dtype=float32)
        dict_bin_meCHHW[chr] = np.zeros(shape, dtype=float32)
        dict_bin_meCHHC[chr] = np.zeros(shape, dtype=float32)
        dict_bin_dp[chr] = np.zeros((nbin,), dtype=float32)
        dict_bin_dpW[chr] = np.zeros((nbin,), dtype=float32)
        dict_bin_dpC[chr] = np.zeros((nbin,), dtype=float32)
        dict_bin_covnCG[chr] = np.zeros(shape, dtype=int32)
        dict_bin_covnCGW[chr] = np.zeros(shape, dtype=int32)
        dict_bin_covnCGC[chr] = np.zeros(shape, dtype=int32)
        dict_bin_covnCHG[chr] = np.zeros(shape, dtype=int32)
        dict_bin_covnCHGW[chr] = np.zeros(shape, dtype=int32)
        dict_bin_covnCHGC[chr] = np.zeros(shape, dtype=int32)
        dict_bin_covnCHH[chr] = np.zeros(shape, dtype=int32)
        dict_bin_covnCHHW[chr] = np.zeros(shape, dtype=int32)
        dict_bin_covnCHHC[chr] = np.zeros(shape, dtype=int32)
        dict_bin_CGden[chr] = np.zeros((nbin,), dtype=float32)
        dict_bin_CGWden[chr] = np.zeros((nbin,), dtype=float32)
        dict_bin_CGCden[chr] = np.zeros((nbin,), dtype=float32)
        dict_bin_CHGden[chr] = np.zeros((nbin,), dtype=float32)
        dict_bin_CHGWden[chr] = np.zeros((nbin,), dtype=float32)
        dict_bin_CHGCden[chr] = np.zeros((nbin,), dtype=float32)
        dict_bin_CHHden[chr] = np.zeros((nbin,), dtype=float32)
        dict_bin_CHHWden[chr] = np.zeros((nbin,), dtype=float32)
        dict_bin_CHHCden[chr] = np.zeros((nbin,), dtype=float32)


        ###########################################################
        #### init chr-wise data
        ###########################################################
        # sum of meth
        dict_chr_meCG[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_meCGW[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_meCGC[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_meCHG[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_meCHGW[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_meCHGC[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_meCHH[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_meCHHW[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_meCHHC[chr] = np.zeros((DP,), dtype=int64)
        # covnCG
        dict_chr_covnCG[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_covnCGW[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_covnCGC[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_covnCHG[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_covnCHGW[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_covnCHGC[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_covnCHH[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_covnCHHW[chr] = np.zeros((DP,), dtype=int64)
        dict_chr_covnCHHC[chr] = np.zeros((DP,), dtype=int64)
        # mean meth
        dict_chr_mmeCG[chr] = np.zeros((DP,), dtype=float32)
        dict_chr_mmeCGW[chr] = np.zeros((DP,), dtype=float32)
        dict_chr_mmeCGC[chr] = np.zeros((DP,), dtype=float32)
        dict_chr_mmeCHG[chr] = np.zeros((DP,), dtype=float32)
        dict_chr_mmeCHGW[chr] = np.zeros((DP,), dtype=float32)
        dict_chr_mmeCHGC[chr] = np.zeros((DP,), dtype=float32)
        dict_chr_mmeCHH[chr] = np.zeros((DP,), dtype=float32)
        dict_chr_mmeCHHW[chr] = np.zeros((DP,), dtype=float32)
        dict_chr_mmeCHHC[chr] = np.zeros((DP,), dtype=float32)
        
    ###########################################################
    #### summarization
    ###########################################################

    S = float32(10_000)
    value: BinCov = BinCov()
    for key, value in dict_binning.items():
        chr, bin = key
        # if chr not in chrs_valid: continue
        if value.length <= 0: continue

        ##### bin-wise
        # meth CG
        dict_bin_meCG[chr][bin, :] = cumratio(value.meCGW[:DP] + value.meCGC[:DP], value.covnCGW[:DP] + value.covnCGC[:DP], float32)/S
        # meth CG W
        dict_bin_meCGW[chr][bin, :] = cumratio(value.meCGW[:DP], value.covnCGW[:DP], float32)/S
        # meth CG C
        dict_bin_meCGC[chr][bin, :] = cumratio(value.meCGC[:DP], value.covnCGC[:DP], float32)/S   
        # meth CHG
        dict_bin_meCHG[chr][bin, :] = cumratio(value.meCHGW[:DP] + value.meCHGC[:DP], value.covnCHGW[:DP] + value.covnCHGC[:DP], float32)/S
        # meth CHG W
        dict_bin_meCHGW[chr][bin, :] = cumratio(value.meCHGW[:DP], value.covnCHGW[:DP], float32)/S
        # meth CHG C
        dict_bin_meCHGC[chr][bin, :] = cumratio(value.meCHGC[:DP], value.covnCHGC[:DP], float32)/S 
        # meth CHH
        dict_bin_meCHH[chr][bin, :] = cumratio(value.meCHHW[:DP] + value.meCHHC[:DP], value.covnCHHW[:DP] + value.covnCHHC[:DP], float32)/S
        # meth CHH W
        dict_bin_meCHHW[chr][bin, :] = cumratio(value.meCHHW[:DP], value.covnCHHW[:DP], float32)/S 
        # meth CHH C
        dict_bin_meCHHC[chr][bin, :] = cumratio(value.meCHHC[:DP], value.covnCHHC[:DP], float32)/S 
        
        # covn of C
        dict_bin_covnCG[chr][bin,:] = cumsumrev(value.covnCGW[:DP] + value.covnCGC[:DP])
        dict_bin_covnCGW[chr][bin,:] = cumsumrev(value.covnCGW[:DP])
        dict_bin_covnCGC[chr][bin,:] = cumsumrev(value.covnCGC[:DP])
        dict_bin_covnCHG[chr][bin,:] = cumsumrev(value.covnCHGW[:DP] + value.covnCHGC[:DP])
        dict_bin_covnCHGW[chr][bin,:] = cumsumrev(value.covnCHGW[:DP])
        dict_bin_covnCHGC[chr][bin,:] = cumsumrev(value.covnCHGC[:DP])
        dict_bin_covnCHH[chr][bin,:] = cumsumrev(value.covnCHHW[:DP] + value.covnCHHC[:DP])
        dict_bin_covnCHHW[chr][bin,:] = cumsumrev(value.covnCHHW[:DP])
        dict_bin_covnCHHC[chr][bin,:] =cumsumrev(value.covnCHHC[:DP])

        # C density
        dict_bin_CGden[chr][bin] = (value.nCGW + value.nCGC)/value.length
        dict_bin_CGWden[chr][bin] = (value.nCGW)/value.length
        dict_bin_CGCden[chr][bin] = (value.nCGC)/value.length
        dict_bin_CHGden[chr][bin] = (value.nCHGW + value.nCHGC)/value.length
        dict_bin_CHGWden[chr][bin] = (value.nCHGW)/value.length
        dict_bin_CHGCden[chr][bin] = (value.nCHGC)/value.length
        dict_bin_CHHden[chr][bin] = (value.nCHHW + value.nCHHC)/value.length
        dict_bin_CHHWden[chr][bin] = (value.nCHHW)/value.length
        dict_bin_CHHCden[chr][bin] = (value.nCHHC)/value.length

        # depth
        dict_bin_dp[chr][bin] = (value.dpW + value.dpW)/value.length
        dict_bin_dpW[chr][bin] = value.dpW/value.length
        dict_bin_dpC[chr][bin] = value.dpC/value.length


        ##### chr-wise
        ## chr-wise me
        dict_chr_meCG[chr] += value.meCGW[:DP] + value.meCGC[:DP]
        dict_chr_meCGW[chr] += value.meCGW[:DP]
        dict_chr_meCGC[chr] += value.meCGC[:DP]
        dict_chr_meCHG[chr] += value.meCHGW[:DP] + value.meCHGC[:DP]
        dict_chr_meCHGW[chr] += value.meCHGW[:DP]
        dict_chr_meCHGC[chr] += value.meCHGC[:DP]
        dict_chr_meCHH[chr] += value.meCHHW[:DP] + value.meCHHC[:DP]
        dict_chr_meCHHW[chr] += value.meCHHW[:DP]
        dict_chr_meCHHC[chr] += value.meCHHC[:DP]
        ## chr-wise cov
        dict_chr_covnCG[chr] += value.covnCGW[:DP] + value.covnCGC[:DP]
        dict_chr_covnCGW[chr] += value.covnCGW[:DP]
        dict_chr_covnCGC[chr] += value.covnCGC[:DP]
        dict_chr_covnCHG[chr] += value.covnCHGW[:DP] + value.covnCHGC[:DP]
        dict_chr_covnCHGW[chr] += value.covnCHGW[:DP]
        dict_chr_covnCHGC[chr] += value.covnCHGC[:DP]
        dict_chr_covnCHH[chr] += value.covnCHHW[:DP] + value.covnCHHC[:DP]
        dict_chr_covnCHHW[chr] += value.covnCHHW[:DP]
        dict_chr_covnCHHC[chr] += value.covnCHHC[:DP]
    
    for chr in chrs_valid:
        dict_chr_mmeCG[chr] = cumratio(dict_chr_meCG[chr], dict_chr_covnCG[chr], float32)/S
        dict_chr_mmeCGW[chr] = cumratio(dict_chr_meCGW[chr], dict_chr_covnCGW[chr], float32)/S
        dict_chr_mmeCGC[chr] = cumratio(dict_chr_meCGC[chr], dict_chr_covnCGC[chr], float32)/S
        dict_chr_mmeCHG[chr] = cumratio(dict_chr_meCHG[chr], dict_chr_covnCHG[chr], float32)/S
        dict_chr_mmeCHGW[chr] = cumratio(dict_chr_meCHGW[chr], dict_chr_covnCHGW[chr], float32)/S
        dict_chr_mmeCHGC[chr] = cumratio(dict_chr_meCHGC[chr], dict_chr_covnCHGC[chr], float32)/S
        dict_chr_mmeCHH[chr] = cumratio(dict_chr_meCHH[chr], dict_chr_covnCHH[chr], float32)/S
        dict_chr_mmeCHHW[chr] = cumratio(dict_chr_meCHHW[chr], dict_chr_covnCHHW[chr], float32)/S
        dict_chr_mmeCHHC[chr] = cumratio(dict_chr_meCHHC[chr], dict_chr_covnCHHC[chr], float32)/S


    ###########################################################
    #### exporting
    ###########################################################
    
    # dict_chr_me = {
    #     'CG':{'double': dict_chr_meCG,'W': dict_chr_meCGW, 'C':dict_chr_meCGC},
    #     'CHG':{'double': dict_chr_meCHG,'W': dict_chr_meCHGW, 'C':dict_chr_meCHGC},
    #     'CHH':{'double': dict_chr_meCHH,'W': dict_chr_meCHHW, 'C':dict_chr_meCHHC}
    #     }
    dict_chr_covn = {
        'CG':{'double': dict_chr_covnCG,'W': dict_chr_covnCGW, 'C':dict_chr_covnCGC},
        'CHG':{'double': dict_chr_covnCHG,'W': dict_chr_covnCHGW, 'C':dict_chr_covnCHGC},
        'CHH':{'double': dict_chr_covnCHH,'W': dict_chr_covnCHHW, 'C':dict_chr_covnCHHC}
        }
    dict_chr_me = {
        'CG':{'double': dict_chr_mmeCG,'W': dict_chr_mmeCGW, 'C':dict_chr_mmeCGC},
        'CHG':{'double': dict_chr_mmeCHG,'W': dict_chr_mmeCHGW, 'C':dict_chr_mmeCHGC},
        'CHH':{'double': dict_chr_mmeCHH,'W': dict_chr_mmeCHHW, 'C':dict_chr_mmeCHHC}
        }
    dict_bin_me = {
        'CG':{'double': dict_bin_meCG,'W': dict_bin_meCGW, 'C':dict_bin_meCGC},
        'CHG':{'double': dict_bin_meCHG,'W': dict_bin_meCHGW, 'C':dict_bin_meCHGC},
        'CHH':{'double': dict_bin_meCHH,'W': dict_bin_meCHHW, 'C':dict_bin_meCHHC}
    }
    dict_bin_covn = {
        'CG':{'double': dict_bin_covnCG,'W': dict_bin_covnCGW, 'C':dict_bin_covnCGC},
        'CHG':{'double': dict_bin_covnCHG,'W': dict_bin_covnCHGW, 'C':dict_bin_covnCHGC},
        'CHH':{'double': dict_bin_covnCHH,'W': dict_bin_covnCHHW, 'C':dict_bin_covnCHHC}
    }
    dict_bin_cden = {
        'CG':{'double': dict_bin_CGden,'W': dict_bin_CGWden, 'C':dict_bin_CGCden},
        'CHG':{'double': dict_bin_CHGden,'W': dict_bin_CHGWden, 'C':dict_bin_CHGCden},
        'CHH':{'double': dict_bin_CHHden,'W': dict_bin_CHHWden, 'C':dict_bin_CHHCden}
    }
    dict_bin_depth = {'double': dict_bin_dp,'W': dict_bin_dpW, 'C':dict_bin_dpC}

    params['dict_chr_me'] = dict_chr_me
    params['dict_chr_covn'] = dict_chr_covn
    params['dict_bin_me'] = dict_bin_me
    params['dict_bin_covn'] = dict_bin_covn
    params['dict_bin_cden'] = dict_bin_cden
    params['dict_bin_depth'] = dict_bin_depth
    return None

def plot_chr_wise_me() -> None:
    chrs = params['chrs_valid']
    dict_chr_me = params['dict_chr_me']
    dict_chr_covn = params['dict_chr_covn']
    DP_valid = params['MAX_DP_BY_FIG'] + 1
    img_dir = params['img_dir']

    strand = 'double'
    for cg in CONTEXTS:
        figwidth = min(12, len(chrs)*0.18 + 0.3)
        ylim = 0
        me_list = []
        dps_list = []
        for dp in range(1, DP_valid):
            # me = [dict_chr_me[cg][strand][chr][dp]/dict_chr_cov[cg][strand][chr][dp] for chr in chrs]
            me = [ dict_chr_me[cg][strand][chr][dp] for chr in chrs ]
            dps = [ dict_chr_covn[cg][strand][chr][dp] for chr in chrs ]
            me_list.append(me)
            dps_list.append(dps)
            ylim = max(ylim, max(me))
        ylim = axlimit(ylim)
        for dp in range(DP_valid-1):
            # [0, DP-1] -> [1, DP+1]
            fig, ax = plt.subplots(figsize = (figwidth, 2))
            dps = np.array(dps_list[dp])
            mes = np.array(me_list[dp])
            i = dps > 0 # remove uncovered chrs
            ax.bar(np.array(chrs)[i], mes[i], color=COLS[0])
            ax.set_ylim(0, ylim)
            plt.xticks(rotation=45)
            ax.set_ylabel('mean DNAme level')
            filename = f'{img_dir}/meth-chr-{cg}-{strand}-dp{dp+1}'
            plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
            if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
            plt.close()
    return None

def plot_binning_meth() -> None:
    chrs = params['chrs_valid']
    dict_bin_me = params['dict_bin_me']
    dict_bin_covn = params['dict_bin_covn']
    binSize = params['binSize']
    img_dir = params['img_dir']
    save_svg = params['save_svg']

    dp = 1 # only plot depth >= 1
    nchr = len(chrs)

    for cg in CONTEXTS:
        for strand in ['double', 'single']:
            if strand == 'double':
                me = dict_bin_me[cg][strand]
                # nan mask for little covered bins
                covn = dict_bin_covn[cg][strand]
                maxbins = max([me[chr].shape[0] for chr in chrs])
                add = int(cg!='CG')*0.2 + int(strand=='single')*0.1
                figheight = nchr*0.35 + 1
                figwidth = maxbins/250 + 1 + add
            else:
                me = dict_bin_me[cg]['W']
                me2 = dict_bin_me[cg]['C']
                covn = dict_bin_covn[cg]['W']
                covn2 = dict_bin_covn[cg]['C']
                maxbins = max([me[chr].shape[0] for chr in chrs])
                figheight = nchr*0.6 + 1
                figwidth = maxbins/250 + 1.2 + add

            fig, axs = plt.subplots(nchr, 1, sharex=True, sharey=True, figsize=(figwidth, figheight))
            if nchr == 1: axs = [axs]
            fig.subplots_adjust(hspace=0)
            
            plt.xlim(-2, prefixBpSize(np.array([maxbins * binSize]))[0]+2)
            kylabel = int((nchr+1) / 2) - 1 # which subplot to place ylabel
            for i, chr in enumerate(chrs):
                shape = np.shape(me[chr])
                x = np.arange(shape[0]) * binSize
                x, prefix = prefixBpSize(x)
                
                # trim outliers in non-CG plots
                y = me[chr][:,dp]
                if strand == 'single':
                    y2 = me2[chr][:,dp]
                if cg != 'CG':
                    y = np.fmin(np.nanquantile(y, 0.99), y)
                    if strand == 'single':
                        y2 = np.fmin(np.nanquantile(y, 0.99), y2)

                # print(np.sum(covn[chr][:,dp] >= 3))
                mask1 = covn[chr][:,dp] < 3
                y[mask1] = nan
                if strand == 'double':
                    axs[i].scatter(x, y, s=1, c=COLS_AREA[0])
                else:
                    y2[covn2[chr][:,dp] < 3] = nan
                    axs[i].scatter(x, y, s=1, c=COLS_AREA[1])
                    axs[i].scatter(x, -y2, s=1, c=COLS_AREA[0])
                    axs[i].hlines(0, 0, max(x), color='#666666', linestyles='dashed', linewidths=1)
                # if i % 2 == nchr % 2:
                if i % 2 != kylabel % 2:
                    axs[i].yaxis.tick_right()
                if i == kylabel: # ensure ylabel on the left always
                    axs[i].yaxis.set_label_position('left')
                    axs[i].set_ylabel('mean DNAme level')
                axs[i].text(axs[i].get_xlim()[1], axs[i].get_ylim()[0], chr, horizontalalignment='right', verticalalignment='bottom')
            # if prefix == '':
            #     axs[i].set_xlabel('genome coordinate')
            # else:
            plt.xlabel(f'genome coordinate ({prefix})')
            # if strand == 'single': axs[i].legend()

            filename = f'{img_dir}/meth-bin-{cg}-{strand}-dp{dp}'
            plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
            if save_svg:
                plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
            plt.close()
    return None

def plot_binning_depth() -> None:
    chrs = params['chrs_valid']
    dict_bin_depth = params['dict_bin_depth']
    binSize = params['binSize']
    img_dir = params['img_dir']
    nchr = len(chrs)
    save_svg = params['save_svg']

    for strand in ['double', 'single']:
        if strand == 'double':
            dep = dict_bin_depth[strand]
            maxbins = max([dep[chr].shape[0] for chr in chrs])
            add = int(strand=='single')*0.1
            figheight = nchr*0.4 + 1
            figwidth = maxbins/200 + 1 + add
        else:
            dep = dict_bin_depth['W']
            dep2 = dict_bin_depth['C']
            maxbins = max([dep[chr].shape[0] for chr in chrs])
            figheight = nchr*0.7 + 1
            figwidth = maxbins/200 + 1.3 + add

        fig, axs = plt.subplots(nchr, 1, sharex=True, sharey=True, figsize=(figwidth, figheight))
        if nchr == 1: axs = [axs]
        fig.subplots_adjust(hspace=0)

        # ylim
        if strand == 'double':
            maxdp = np.quantile(np.hstack([dep[chr] for chr in chrs]), 0.98)
        else:
            maxdp = np.quantile(np.hstack([np.hstack([dep[chr], dep2[chr]]) for chr in chrs]), 0.98)
            
        kylabel = int((nchr+1) / 2) - 1 # which subplot to place ylabel
        for i, chr in enumerate(chrs):
            shape = np.shape(dep[chr])
            x0 = np.arange(shape[0]) * binSize
            x, prefix = prefixBpSize(x0)
            
            # trim outliers
            y = dep[chr]
            y = np.fmin(maxdp, y)
            if strand == 'single':
                y2 = dep2[chr]
                y2 = np.fmin(maxdp, y2)
            
            if strand == 'double':
                axs[i].scatter(x, y, s=1, c=COLS_AREA[0])
            else:
                axs[i].scatter(x, y, s=1, c=COLS_AREA[1])
                axs[i].scatter(x, -y2, s=1, c=COLS_AREA[0])
                axs[i].hlines(0, 0, max(x), color='#666666', linestyles='dashed', linewidths=1)
            # if i%2 == nchr%2:
            if i % 2 != kylabel % 2:
                axs[i].yaxis.tick_right()
            if i == kylabel: # ensure ylabel on the left always
                axs[i].yaxis.set_label_position('left')
                axs[i].set_ylabel('mean read depth')
            axs[i].text(axs[i].get_xlim()[1], axs[i].get_ylim()[0], chr, horizontalalignment='right', verticalalignment='bottom')
        # if prefix == '':
        #     axs[i].set_xlabel('genome coordinate')
        # else:
        plt.xlabel(f'genome coordinate ({prefix})')

        filename = f'{img_dir}/depth-bin-{strand}'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        if save_svg:
            plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    return None
    