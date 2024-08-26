
import math
import os
from typing import NamedTuple, List, Tuple, Dict, Any

from src.utils import myMakeDirs
from src.utils import getBins


# class Parameters():
#     fafile: str
#     fa: MyFastaFile
#     bamfile: str
#     # out_atcg: str
#     # out_cg: str
#     # out_bed: str
#     # chr: str
#     start: int
#     end: int
#     step: int
#     quality_threshold: int
#     context_size: int
#     cgkmer: int
#     coordinate_base: int  # 0/1-based
#     swap_strand: bool # swap the counts of forward and reverse strands
#     read_quality: int

#     chr_MT: str
#     chr_lambda: str
#     chr_plastid: str
#     testchrs: str
#     binSize:int # bin length
#     bins_lambda: int # #bins
#     bins_MT: int
#     MAXDEPTH: int # max depth for counting
    
#     # gene pading for pan-gene summary
#     gene_padding: int # up-/down-stream
#     gene_breaks: int # breaks of up body down

#     # read sampling
#     reads_to_sample: int
#     reads_sample_window: int # bp
#     max_reads_per_batch: int

#     # pangene meth
#     PANGENE_SAMPLED: int

#     binSize_lambda: Tuple
#     binSize_MT: Tuple
#     binSizeContig: Dict
#     binsContig: Dict
#     chrs_alias: Dict




##############################################

# def config_params(args: NamedTuple):
def config_params():
    
    fafile = 'genome/hg38.fa'
    # fafile = 'genome/GCF_000001735.4_TAIR10.1_genomic.fna'
    bamfile = 'data/ENCFF873NOV.bam'
    # bamfile = "data/SRR19264609.sorted.bam"
    # bamfile = 'data/ENCFF000LTD.sorted.bam'
    # bamfile = 'data/SRR23252523.sorted.bam'

    # chr_MT = 'NC_037304.1'
    # chr_lambda = 'NC_000932.1'
    chr_MT = 'chrM'
    chr_lambda = 'chrL'
    # chr_plastid = 'NC_000932.1'
    chr_plastid = ''

    gtffile = "genome/gencode.v46.annotation.gtf.gz"
    # gtffile = "/home/yance/bsDoctor/genome/GCF_000001735.4_TAIR10.1_genomic.gtf.gz"

    # testchrs = 'all'
    # testchrs = 'chr21,chr22'
    testchrs = 'chr1,chr2,chr3,chr18,chr19,chr20,chr21,chr22,chrX'
    # testchrs = 'NC_003070.9,NC_003071.7,NC_003074.8,NC_003075.7,NC_003076.8'
    # testchrs = 'NC_003075.7,NC_003076.8'

    binSize = 100_000 # bin length
    bins_lambda = 1000 # #bins
    bins_MT = 1000

    # max depth for counting
    # depth = max{depth, MAXDEPTH}
    MAXDEPTH = 200

    # gene pading for pan-gene summary
    gene_padding = 5_000 # up-/down-stream
    gene_breaks = 50 # breaks of up body down

    # quality read sampling
    reads_to_sample = 100_00
    # reads_to_sample = 100_000
    reads_sample_window = 10 # bp
    max_reads_per_batch = 500

    # pangene meth
    PANGENE_SAMPLED = 100


    params = dict()
    
    params['fafile'] = fafile
    params['bamfile'] = bamfile
    params['gtffile'] = gtffile

    params['start'] = 0
    params['end'] = math.inf
    params['step'] = 2000
    params['quality_threshold'] = 0
    params['read_quality'] = 0
    params['nuclear_sampling_step'] = 100
    params['nuclear_sampling_spacing'] = 10_000
    params['context_size'] = 3
    params['cgkmerSize'] = 6
    params['coordinate_base'] = 1
    params['swap_strand'] = False
    params['chr_MT'] = chr_MT
    params['chr_lambda'] = chr_lambda
    params['chr_plastid'] = chr_plastid
    params['testchrs'] = testchrs
    params['binSize'] = binSize
    params['bins_lambda'] = bins_lambda
    params['bins_MT'] = bins_MT
    params['MAXDEPTH'] = MAXDEPTH
    # params['meth_bins'] = 20 #  #bins of meth
    params['gene_padding'] = gene_padding
    params['gene_breaks'] = gene_breaks
    params['reads_to_sample'] = reads_to_sample
    params['reads_sample_window'] = reads_sample_window
    params['max_reads_per_batch'] = max_reads_per_batch
    params['PANGENE_SAMPLED'] = PANGENE_SAMPLED
    params['MAX_COORDINATE'] = 2**63-1
    params['MAX_DP_BY_FIG'] = 30 # max dp for which figures plotted by each dp
    params['MAX_DP_CG_MOTIF'] = 50
    params['work_dir'] = '.'
    # params['report_dir'] = os.path.join(params['work_dir'], 'report')
    params['report_dir'] = params['work_dir']
    params['img_dir'] = os.path.join(params['report_dir'], 'img2')
    
    myMakeDirs(params['work_dir'])
    myMakeDirs(params['report_dir'])
    myMakeDirs(params['img_dir'])

    return params


global params
params = config_params()

    # params = Parameters(
    #     fafile=fafile, 
    #     fa=None,
    #     bamfile=bamfile,
    #     # out_atcg='', # three types of outputs
    #     # out_cg='',
    #     # out_bed='',
    #     # chr='chr20',  # 'all' for all chrs
    #     start=0,
    #     end=math.inf,
    #     # end=50_000_000,
    #     step=2000,
    #     quality_threshold=0, # base seq quality
    #     read_quality=0, # read mapping quality threshold
    #     context_size=3, # size of CHG/...
    #     cgkmer=6,
    #     coordinate_base=1, # 0/1-based
    #     swap_strand=False,  # swap counts of two strands

    #     chr_MT = chr_MT,
    #     chr_lambda = chr_lambda,
    #     chr_plastid = chr_plastid,
    #     testchrs = testchrs,
    #     binSize = binSize,
    #     bins_lambda = bins_lambda,
    #     bins_MT = bins_MT,
    #     MAXDEPTH = MAXDEPTH,
    #     gene_padding = gene_padding,
    #     gene_breaks = gene_breaks,
    #     reads_to_sample = reads_to_sample,
    #     reads_sample_window = reads_sample_window,
    #     max_reads_per_batch = max_reads_per_batch,
    #     PANGENE_SAMPLED = PANGENE_SAMPLED,
    #     binSize_lambda = None,
    #     binSize_MT = None,
    #     binSizeContig = None,
    #     binsContig = None,
    #     chrs_alias = None
    #     )


# data exported to Jinja
global data
data = {}


