#! /bin/env python

import re
import pickle
# import numpy as np

from jinja2 import Environment,FileSystemLoader

#######

from src.config import *
from src.utils import *

from src.coverage import *
from src.quality import *
from src.pangene import *
from src.summary import *
from src.updateBinning import *
from src.mt import *
from src.lambdaDNA import *
from src.plastid import *
from src.whole_genome import *
from src.bin_chr_wise import *
from src.CpG_strandness import *
from src.me_bias import *
from src.depth_bias import *
from src.CpG_motif import *
from src.RRBS import *
from src.saturation_curve_and_DNA_lost import *

#############################
#### config parameters
#############################

def config_params_further() -> None:

    SAMPLE = re.search(r'/?([^/\.]+)\.?[^/]*$', params['bamfile']).groups()[0]
    # print(SAMPLE)
    data['SAMPLE'] = SAMPLE
    params['SAMPLE'] = SAMPLE
    data['title'] = f'bsDoctor Report of Sample: {SAMPLE}'

    fa = MyFastaFile(params['fafile'])
    params['fa'] = fa

    reference_length = fa.get_reference_length
    params['reference_length'] = reference_length

    bam = MyAlignmentFile(params['bamfile'], 'rb')
    params['bam'] = bam

    # chrs to sample
    # chrs_valid = chrs.split(',')
    chrs = params['testchrs']
    if chrs == 'all':
        # chrs_valid = [
        #     chr for chr in fa.references
        #     if not chr.endswith(('_random', '_alt')) and not chr.startswith('chrUn_') and chr not in (params['chr_lambda'], params['chr_MT'], params['chr_plastid'])
        # ]
        # order chrs according to bam file
        chrs_valid = [c.contig for c in bam.get_index_statistics()]
        chrs_valid = excludeContigs(chrs_valid, params['contigs_ex_start'], params['contigs_ex_end'])
        chrs_valid = [
            c for c in chrs_valid if c in fa.references and c not in (params['chr_lambda'], params['chr_MT'], params['chr_plastid'])
        ]
    else:
        chrs_valid = chrs.split(',')

    # print(chrs_valid)

    # lengths of chrs
    lens_valid = [reference_length(chr) for chr in chrs_valid]
    params['chrs_valid'] = chrs_valid
    params['lens_valid'] = lens_valid
    # params['max_chr_len'] = max(lens_valid)

    # nuclear bin size
    if params['binSize'] == 0:
        binSize = getBins(max(lens_valid), params['nbins_chr'])[0]
        params['binSize'] = binSize
    
    # sampling size
    binSize = params['binSize']
    prop = params['sampling_prop']
    if params['binSize'] <= 5000 or prop > 1-1e-6:
        params['nuclear_sampling_step'] = binSize
        params['nuclear_sampling_spacing'] = 0
    else:
        samples_in_a_bin = 10
        step_size = max(1000, int(binSize*prop/samples_in_a_bin))
        params['nuclear_sampling_step'] = step_size
        params['nuclear_sampling_spacing'] = min(int(step_size*(1/prop-1)), int(binSize*(1-prop)/samples_in_a_bin))
        if params['nuclear_sampling_spacing'] < step_size/3:
            params['nuclear_sampling_step'] = step_size
            params['nuclear_sampling_spacing'] = 0

    print(f'nuclear bin size: {params["binSize"]}; sampling: {params["nuclear_sampling_step"]}, {params["nuclear_sampling_spacing"]}')

    # prefix
    a, b = prefixBpSize(params['binSize'])
    if np.abs(a - np.ceil(a)) < 1e-6:
        data['bin_size'] = f'{int(a)} {b}'
    else:
        data['bin_size'] = f'{a} {b}'
    # data['bin_size']

    chr_MT = params['chr_MT']
    chr_lambda = params['chr_lambda']
    chr_plastid = params['chr_plastid']
    chrs_alias = {'lambda': chr_lambda, 'MT': chr_MT, 'plastid': chr_plastid}
    params['chrs_alias'] = chrs_alias

    bins_lambda = params['bins_lambda']
    bins_mt = params['bins_lambda']
    bins_plastid = params['bins_plastid']

    binSizeContig ={}
    binsContig = {}
    if data['include_lambda']:
        binSize_lambda = getBins(reference_length(chr_lambda), bins_lambda)
        binSizeContig[chr_lambda] = binSize_lambda[0]
        binsContig[chr_lambda] = binSize_lambda[1]
        params['binSize_lambda'] = binSize_lambda
        data['binSize_lambda'] = fi(binSize_lambda[0])
    if data['include_mt']:
        binSize_MT = getBins(reference_length(chr_MT), bins_mt)
        binSizeContig[chr_MT] = binSize_MT[0]
        binsContig[chr_MT] = binSize_MT[1]
        params['binSize_MT'] = binSize_MT
        data['binSize_MT'] = fi(binSize_MT[0])
    if data['include_plastid']:
        binSize_plastid = getBins(reference_length(chr_plastid), bins_plastid)
        binSizeContig[chr_plastid] = binSize_plastid[0]
        binsContig[chr_plastid] = binSize_plastid[1]
        params['binSize_plastid'] = binSize_plastid
        data['binSize_plastid'] = fi(binSize_plastid[0])

    params['binSizeContig'] = binSizeContig
    params['binsContig'] = binsContig
    
    print(f'{binSizeContig=}', f'{binsContig=}')

    # ('chr1', 23) -> BinCov()
    dict_binning = dict()
    params['dict_binning'] = dict_binning

    # CpG-motif
    if data['include_motif']:
        # 'NNCGNN' -> KmerCov()
        dict_cgkmer = dict()
        params['dict_cgkmer'] = dict_cgkmer

    params['reads_to_sample'] = min([params['reads_to_sample'], bam.mapped])

    return None
    

#######################################
#### summarization and plots
#######################################

def compute_and_plot():

    # DNA content
    compt_plot_DNA_content()

    # base and read quality
    if data['include_quality']:
        # quality stats       
        quality = Quality()
        params['quality'] = quality
        compt_quality()
        plot_read_length()
        plot_base_quality()
        plot_read_map_quality()
        plot_bar_base_cigar()

    # pangene
    if data['include_pangene']:
        pangene_sampling()
        pangene_compt_plot_meth()

    # MT
    if data['include_mt']:
        compt_MT()
        if data['mt_is_covered'] == 1:
            plot_mt_depth_binning()
            plot_mt_base_error_rate()

    # lambda
    if data['include_lambda']:
        compt_lambda()
        if data['lambda_is_covered'] == 1:
            plot_lambda_depth_binning()
            plot_lambda_base_error_rate()

    # plastid
    if data['include_plastid']:
        compt_plastid()
        if data['plastid_is_covered'] == 1:
            plot_plastid_depth_binning()
            plot_plastid_base_error_rate()

    # nuclear binning
    nuclear_sampling()
    compt_chr_and_bin_wise()
    plot_chr_wise_me()

    # nuclear whole-genome
    compt_whole_genome()
    plot_base_error_rate_by_AT()
    plot_theroretical_me_bias()
    plot_hist_me()

    # CpG strandness
    plot_heatmap_stranded_CpG_depth()
    plot_bar_CpG_depth_difference()
    plot_bar_double_srtanded_cpg()
    plot_heatmap_stranded_CpG_meth()
    plot_heatmap_stranded_meth_diff()

    # binning
    plot_binning_meth()
    plot_binning_depth()

    # me bias
    plot_me_vs_depth_eq_k()
    plot_me_vs_depth_ge_k()
    plot_me_vs_missing()
    plot_me_and_covrate_vs_cytosine_density()

    # depth bias
    plot_depth_vs_cytosine_density()
    plot_depth_dist_by_low_high_me()
    plot_covrate_vs_depth_of_cytosine()
    plot_depth_dist_of_cytosine()
    plot_covrate_vs_depth_of_whole_genome()
    plot_depth_dist_of_whole_genome()
    plot_depth_watson_vs_crick()
    # plot_depth_overall_vs_me()
    plot_depth_AT_vs_CG()

    # CpG-motif
    if data['include_motif']:
        compt_CpG_motif()
        plot_hist_CpG_motif_freq()
        plot_CpG_motif_freq_watson_vs_crick()
        plot_CpG_motif_me_watson_vs_crick()
        plot_CpG_motif_depth_vs_freq()
        plot_CpG_motif_me_vs_freq()
        plot_CpG_motif_covrate_vs_depth()
        plot_CpG_motif_me_vs_depth()
        # RRBS
        compt_RRBS()
        plot_RRBS_CpG_motifs()

    # saturation curve and DNA lost
    if data['include_saturation']:
        compt_plot_DNA_composition()
        plot_saturation_curve()

#######################################
#### write report
#######################################

def write_report() -> None:
    # with open('rrbs-human2/data.pickle', 'rb') as f:
    #     data = pickle.load(f)
    
    env = Environment(loader=FileSystemLoader('report/'))
    env.globals["include_file"] = include_file
    template = env.get_template('main.html')
    # print(data)
    temp_out = template.render(alldata=data)   
    html_out = f"{data['report_dir']}/report-{data['SAMPLE']}.html"

    with open(html_out, 'w', encoding='utf-8') as f:
        f.writelines(temp_out)
    return None

if __name__ == "__main__":
    desc = "bsDoctor: Quality Diagnosis for Bisulfite-Seq Data"
    parser = MyArgumentParser(description=desc)
    options = parser.parse_args()
    
    check_args(options)
    config_params(options)      
    config_params_further()

    print(data)
    compute_and_plot()

    with open(f"{params['report_dir']}/data.pickle", 'wb') as fd:
        pickle.dump(data, fd)

    write_report()
        