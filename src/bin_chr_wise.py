
import tqdm
from src.utils import *
from src.config import params, data
from src.coverage import CovLambda, GenomicIntervalGenerator
from src.updateBinning import update_cgkmer, update_binning


def nuclear_sampling():
    fa = params['fa']
    bam = params['bam']
    testchrs = params['testchrs']
    step = params['nuclear_sampling_step']
    spacing = params['nuclear_sampling_spacing']

    intervals = iter(GenomicIntervalGenerator(
        fa, 
        chrs= testchrs,
        # chrs='all',
        start = 0,    
        end = params['MAX_COORDINATE'],
        step=step, # 1_000
        spacing=spacing #10_000
        ))

    intervals_list = list(intervals)
    for i in tqdm.trange(len(intervals_list)):
        detailedIntvl = bam.detailedCoverage(intervals_list[i])
        update_cgkmer(detailedIntvl)
        update_binning(detailedIntvl)

def compt_chr_and_bin_wise():
    fa = params['fa']
    reference_length = params['reference_length']
    binSize = params['binSize']
    dict_binning = params['dict_binning']
    chrs_valid = params['chrs_valid']
    DP_valid = params['MAX_DP_BY_FIG']

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
    dict_bin_dp = dict()
    dict_bin_dpW = dict()
    dict_bin_dpC = dict()
    # DP_valid = 30
    

    for i, chr in enumerate(chrs_valid):
        nbin = int(np.ceil(reference_length(chr)/binSize))
        dict_bin_meCG[chr] = np.zeros((nbin, DP_valid), dtype=float)
        dict_bin_meCGW[chr] = np.zeros((nbin, DP_valid), dtype=float)
        dict_bin_meCGC[chr] = np.zeros((nbin, DP_valid), dtype=float)
        dict_bin_meCHG[chr] = np.zeros((nbin, DP_valid), dtype=float)
        dict_bin_meCHGW[chr] = np.zeros((nbin, DP_valid), dtype=float)
        dict_bin_meCHGC[chr] = np.zeros((nbin, DP_valid), dtype=float)
        dict_bin_meCHH[chr] = np.zeros((nbin, DP_valid), dtype=float)
        dict_bin_meCHHW[chr] = np.zeros((nbin, DP_valid), dtype=float)
        dict_bin_meCHHC[chr] = np.zeros((nbin, DP_valid), dtype=float)
        dict_bin_dp[chr] = np.zeros((nbin, ), dtype=float)
        dict_bin_dpW[chr] = np.zeros((nbin, ), dtype=float)
        dict_bin_dpC[chr] = np.zeros((nbin, ), dtype=float)

    ###########################################################
    #### init chr-wise data
    ###########################################################

    # sum of meth
    dict_chr_meCG = dict.fromkeys(chrs_valid, 0)
    dict_chr_meCGW = dict.fromkeys(chrs_valid, 0)
    dict_chr_meCGC = dict.fromkeys(chrs_valid, 0)
    dict_chr_meCHG = dict.fromkeys(chrs_valid, 0)
    dict_chr_meCHGW = dict.fromkeys(chrs_valid, 0)
    dict_chr_meCHGC = dict.fromkeys(chrs_valid, 0)
    dict_chr_meCHH = dict.fromkeys(chrs_valid, 0)
    dict_chr_meCHHW = dict.fromkeys(chrs_valid, 0)
    dict_chr_meCHHC = dict.fromkeys(chrs_valid, 0)
    dict_chr_meCHG = dict.fromkeys(chrs_valid, 0)
    dict_chr_meCHGW = dict.fromkeys(chrs_valid, 0)
    dict_chr_meCHGC = dict.fromkeys(chrs_valid, 0)
    # covnCG
    dict_chr_covnCG = dict.fromkeys(chrs_valid, 0)
    dict_chr_covnCGW = dict.fromkeys(chrs_valid, 0)
    dict_chr_covnCGC = dict.fromkeys(chrs_valid, 0)
    dict_chr_covnCHG = dict.fromkeys(chrs_valid, 0)
    dict_chr_covnCHGW = dict.fromkeys(chrs_valid, 0)
    dict_chr_covnCHGC = dict.fromkeys(chrs_valid, 0)
    dict_chr_covnCHH = dict.fromkeys(chrs_valid, 0)
    dict_chr_covnCHHW = dict.fromkeys(chrs_valid, 0)
    dict_chr_covnCHHC = dict.fromkeys(chrs_valid, 0)
    dict_chr_covnCHG = dict.fromkeys(chrs_valid, 0)
    dict_chr_covnCHGW = dict.fromkeys(chrs_valid, 0)
    dict_chr_covnCHGC = dict.fromkeys(chrs_valid, 0)

    ###########################################################
    #### traverse sampled intervals
    ###########################################################
    for key, value in dict_binning.items():
        chr, bin = key
        if chr not in chrs_valid: continue

        ##### bin-wise
        # meth CG
        index = value.covnCG[:DP_valid] >= 1
        dict_bin_meCG[chr][bin, index] = value.meCG[:DP_valid][index]/value.covnCG[:DP_valid][index]
        dict_bin_meCG[chr][bin, ~index] = np.nan
        # meth CG W
        index = value.covnCGW[:DP_valid] >= 1
        dict_bin_meCGW[chr][bin, index] = value.meCGW[:DP_valid][index]/value.covnCGW[:DP_valid][index]
        dict_bin_meCGW[chr][bin, ~index] = np.nan
        # meth CG C
        index = value.covnCGC[:DP_valid] >= 1
        dict_bin_meCGC[chr][bin, index] = value.meCGC[:DP_valid][index]/value.covnCGC[:DP_valid][index]
        dict_bin_meCGC[chr][bin, ~index] = np.nan
        # meth CHG
        index = value.covnCHG[:DP_valid] >= 1
        dict_bin_meCHG[chr][bin, index] = value.meCHG[:DP_valid][index]/value.covnCHG[:DP_valid][index]
        dict_bin_meCHG[chr][bin, ~index] = np.nan
        # meth CHG W
        index = value.covnCHGW[:DP_valid] >= 1
        dict_bin_meCHGW[chr][bin, index] = value.meCHGW[:DP_valid][index]/value.covnCHGW[:DP_valid][index]
        dict_bin_meCHGW[chr][bin, ~index] = np.nan
        # meth CHG C
        index = value.covnCHGC[:DP_valid] >= 1
        dict_bin_meCHGC[chr][bin, index] = value.meCHGC[:DP_valid][index]/value.covnCHGC[:DP_valid][index]
        dict_bin_meCHGC[chr][bin, ~index] = np.nan
        # meth CHH
        index = value.covnCHH[:DP_valid] >= 1
        dict_bin_meCHH[chr][bin, index] = value.meCHH[:DP_valid][index]/value.covnCHH[:DP_valid][index]
        dict_bin_meCHH[chr][bin, ~index] = np.nan
        # meth CHH W
        index = value.covnCHHW[:DP_valid] >= 1
        dict_bin_meCHHW[chr][bin, index] = value.meCHHW[:DP_valid][index]/value.covnCHHW[:DP_valid][index]
        dict_bin_meCHHW[chr][bin, ~index] = np.nan
        # meth CHH C
        index = value.covnCHHC[:DP_valid] >= 1
        dict_bin_meCHHC[chr][bin, index] = value.meCHHC[:DP_valid][index]/value.covnCHHC[:DP_valid][index]
        dict_bin_meCHHC[chr][bin, ~index] = np.nan
        # depth
        if value.length < 0: continue
        dict_bin_dp[chr][bin] = value.dp/value.length
        dict_bin_dpW[chr][bin] = value.dpW/value.length
        dict_bin_dpC[chr][bin] = value.dpC/value.length

        ##### chr-wise
        ## chr-wise me
        dict_chr_meCG[chr] += value.meCG[:DP_valid]
        dict_chr_meCGW[chr] += value.meCGW[:DP_valid]
        dict_chr_meCGC[chr] += value.meCGC[:DP_valid]
        dict_chr_meCHG[chr] += value.meCHG[:DP_valid]
        dict_chr_meCHGW[chr] += value.meCHGW[:DP_valid]
        dict_chr_meCHGC[chr] += value.meCHGC[:DP_valid]
        dict_chr_meCHH[chr] += value.meCHH[:DP_valid]
        dict_chr_meCHHW[chr] += value.meCHHW[:DP_valid]
        dict_chr_meCHHC[chr] += value.meCHHC[:DP_valid]
        ## chr-wise cov
        dict_chr_covnCG[chr] += value.covnCG[:DP_valid]
        dict_chr_covnCGW[chr] += value.covnCGW[:DP_valid]
        dict_chr_covnCGC[chr] += value.covnCGC[:DP_valid]
        dict_chr_covnCHG[chr] += value.covnCHG[:DP_valid]
        dict_chr_covnCHGW[chr] += value.covnCHGW[:DP_valid]
        dict_chr_covnCHGC[chr] += value.covnCHGC[:DP_valid]
        dict_chr_covnCHH[chr] += value.covnCHH[:DP_valid]
        dict_chr_covnCHHW[chr] += value.covnCHHW[:DP_valid]
        dict_chr_covnCHHC[chr] += value.covnCHHC[:DP_valid]
    
    ###########################################################
    #### exporting
    ###########################################################
    
    dict_chr_me = {
        'CG':{'double': dict_chr_meCG,'W': dict_chr_meCGW, 'C':dict_chr_meCGC},
        'CHG':{'double': dict_chr_meCHG,'W': dict_chr_meCHGW, 'C':dict_chr_meCHGC},
        'CHH':{'double': dict_chr_meCHH,'W': dict_chr_meCHHW, 'C':dict_chr_meCHHC}
        }
    dict_chr_cov = {
        'CG':{'double': dict_chr_covnCG,'W': dict_chr_covnCGW, 'C':dict_chr_covnCGC},
        'CHG':{'double': dict_chr_covnCHG,'W': dict_chr_covnCHGW, 'C':dict_chr_covnCHGC},
        'CHH':{'double': dict_chr_covnCHH,'W': dict_chr_covnCHHW, 'C':dict_chr_covnCHHC}
        }
    dict_bin_me = {
        'CG':{'double': dict_bin_meCG,'W': dict_bin_meCGW, 'C':dict_bin_meCGC},
        'CHG':{'double': dict_bin_meCHG,'W': dict_bin_meCHGW, 'C':dict_bin_meCHGC},
        'CHH':{'double': dict_bin_meCHH,'W': dict_bin_meCHHW, 'C':dict_bin_meCHHC}
    }
    dict_bin_depth = {'double': dict_bin_dp,'W': dict_bin_dpW, 'C':dict_bin_dpC}

    params['dict_chr_me'] = dict_chr_me
    params['dict_chr_cov'] = dict_chr_cov
    params['dict_bin_me'] = dict_bin_me
    params['dict_bin_depth'] = dict_bin_depth

def plot_chr_wise_me():
    chrs_plot = params['chrs_valid']
    dict_chr_me = params['dict_chr_me']
    dict_chr_cov = params['dict_chr_cov']
    DP_valid = params['MAX_DP_BY_FIG']
    img_dir = params['img_dir']

    strand = 'double'
    for cg in CONTEXTS:
        figwidth = min(12, len(chrs_plot)*0.15 + 1)
        ylim = 0
        me_list = []
        for dp in range(DP_valid):
            me = [dict_chr_me[cg][strand][chr][dp]/dict_chr_cov[cg][strand][chr][dp] for chr in chrs_plot]
            me_list.append(me)
            ylim = max(ylim, max(me))
        ylim = axlimit(ylim)
        for dp in range(DP_valid):
            fig, ax = plt.subplots(figsize = (figwidth,2))
            ax.bar(chrs_plot, me_list[dp], color=COLS[0])
            ax.set_ylim(0, ylim)
            plt.xticks(rotation=45)
            ax.set_ylabel('mean methylation')
            filename = f'{img_dir}/meth-chr-{cg}-{strand}-dp{dp+1}'
            plt.savefig(filename+'.png', transparent=True, dpi=120, bbox_inches='tight')
            plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
            plt.close()

def plot_binning_meth():
    chrs_plot = params['chrs_valid']
    dict_bin_me = params['dict_bin_me']
    binSize = params['binSize']
    img_dir = params['img_dir']

    dp = 0 # only plot depth >= 1
    nchr = len(chrs_plot)

    for cg in CONTEXTS:
        for strand in ['double', 'single']:
            if strand == 'double':
                me = dict_bin_me[cg][strand]
                maxbins = max([me[chr].shape[0] for chr in chrs_plot])
                add = int(cg!='CG')*0.2 + int(strand=='single')*0.1
                figheight = nchr*0.35 + 1
                figwidth = maxbins/250 + 1 + add
            else:
                me = dict_bin_me[cg]['W']
                me2 = dict_bin_me[cg]['C']
                maxbins = max([me[chr].shape[0] for chr in chrs_plot])
                figheight = nchr*0.6 + 1
                figwidth = maxbins/250 + 1.2 + add

            fig, axs = plt.subplots(nchr, 1, sharex=True, sharey=True, figsize=(figwidth, figheight))
            fig.subplots_adjust(hspace=0)
            
            plt.xlim(-2, prefixBpSize(np.array([maxbins * binSize]))[0]+2)
            for i, chr in enumerate(chrs_plot):
                if i + 1 == int((nchr+1) / 2):
                    axs[i].set_ylabel('mean DNAme level')
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
                
                if strand == 'double':
                    axs[i].scatter(x, y, s=1, c=COLS_AREA[0])
                else:
                    axs[i].scatter(x, y, s=1, c=COLS_AREA[1])
                    axs[i].scatter(x, -y2, s=1, c=COLS_AREA[0])
                    axs[i].hlines(0, 0, max(x), color='#666666', linestyles='dashed', linewidths=1)
                if i % 2 == nchr % 2: # ensure ylabel on the left always
                    axs[i].yaxis.set_label_position('right')
                    axs[i].yaxis.tick_right()
                axs[i].text(axs[i].get_xlim()[1], axs[i].get_ylim()[0], chr, horizontalalignment='right', verticalalignment='bottom')
            # if prefix == '':
            #     axs[i].set_xlabel('genome coordinate')
            # else:
            plt.xlabel(f'genome coordinate ({prefix})')
            # if strand == 'single': axs[i].legend()

            filename = f'{img_dir}/meth-bin-{cg}-{strand}-dp{dp+1}'
            plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
            plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
            plt.close()

def plot_binning_depth():
    chrs_plot = params['chrs_valid']
    dict_bin_depth = params['dict_bin_depth']
    binSize = params['binSize']
    img_dir = params['img_dir']
    nchr = len(chrs_plot)

    for strand in ['double', 'single']:
        if strand == 'double':
            dep = dict_bin_depth[strand]
            maxbins = max([dep[chr].shape[0] for chr in chrs_plot])
            add = int(strand=='single')*0.1
            figheight = nchr*0.4 + 1
            figwidth = maxbins/200 + 1 + add
        else:
            dep = dict_bin_depth['W']
            dep2 = dict_bin_depth['C']
            maxbins = max([dep[chr].shape[0] for chr in chrs_plot])
            figheight = nchr*0.7 + 1
            figwidth = maxbins/200 + 1.3 + add

        fig, axs = plt.subplots(nchr, 1, sharex=True, sharey=True, figsize=(figwidth, figheight))
        fig.subplots_adjust(hspace=0)

        # ylim
        if strand == 'double':
            maxdp = np.quantile(np.hstack([dep[chr] for chr in chrs_plot]), 0.98)
        else:
            maxdp = np.quantile(np.hstack([np.hstack([dep[chr], dep2[chr]]) for chr in chrs_plot]), 0.98)
        for i, chr in enumerate(chrs_plot):
            if i + 1 == int((nchr+1) / 2):
                axs[i].set_ylabel('mean read depth')
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
            if i%2 == nchr%2: # ensure ylabel on the left always
                axs[i].yaxis.set_label_position('right')
                axs[i].yaxis.tick_right()
            axs[i].text(axs[i].get_xlim()[1], axs[i].get_ylim()[0], chr, horizontalalignment='right', verticalalignment='bottom')
        # if prefix == '':
        #     axs[i].set_xlabel('genome coordinate')
        # else:
        plt.xlabel(f'genome coordinate ({prefix})')

        filename = f'{img_dir}/depth-bin-{strand}'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    