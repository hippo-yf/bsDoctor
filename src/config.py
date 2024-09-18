
import math
import os
from argparse import ArgumentParser, Namespace
from datetime import datetime

from src.utils import *

import matplotlib
matplotlib.rcParams["svg.fonttype"] = 'none'

__version__ = '1.0.0'


# data in python
global params
params = dict()

# data exported to Jinja template
global data
data = {}


##############################################

# def config_params(args: NamedTuple):
def config_params(options: Namespace = Namespace()) -> None:
    
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
    testchrs = 'chr22'
    # testchrs = 'chr1,chr2,chr3,chr18,chr19,chr20,chr21,chr22,chrX'
    # testchrs = 'NC_003070.9,NC_003071.7,NC_003074.8,NC_003075.7,NC_003076.8'
    # testchrs = 'NC_003075.7,NC_003076.8'

    binSize = 100_000 # bin length
    bins_lambda = 1000 # #bins
    bins_MT = 1000
    bins_plastid = 1000

    # max depth for counting
    # depth = max{depth, MAXDEPTH}
    MAXDEPTH = 200

    # gene pading for pan-gene summary
    gene_padding = 5_000 # up-/down-stream
    gene_breaks = 50 # breaks of up body down

    # quality read sampling
    reads_to_sample = 10_000
    # reads_to_sample = 100_000
    reads_sample_window = 10 # bp
    max_reads_per_batch = 500

    # pangene meth
    PANGENE_SAMPLED = 100
    
    params['fafile'] = fafile
    params['bamfile'] = bamfile
    params['gtffile'] = gtffile

    params['start'] = 0
    params['end'] = math.inf
    # params['step'] = 2000
    params['quality_threshold'] = 0
    params['read_quality'] = 0
    params['nuclear_sampling_step'] = 1000
    params['nuclear_sampling_spacing'] = 10_000
    params['context_size'] = 3
    params['cgkmerSize'] = 6
    # params['coordinate_base'] = 1
    params['swap_strand'] = False
    params['chr_MT'] = chr_MT
    params['chr_lambda'] = chr_lambda
    params['chr_plastid'] = chr_plastid
    params['testchrs'] = testchrs
    params['binSize'] = binSize
    params['bins_lambda'] = bins_lambda
    params['bins_MT'] = bins_MT
    params['bins_plastid'] = bins_plastid
    params['MAXDEPTH'] = MAXDEPTH
    # params['meth_bins'] = 20 #  #bins of meth
    
    # pangene
    params['PANGENE_SAMPLED'] = PANGENE_SAMPLED
    params['gene_padding'] = gene_padding
    params['gene_breaks'] = gene_breaks

    # quality
    params['reads_to_sample'] = reads_to_sample
    params['reads_sample_window'] = reads_sample_window
    params['max_reads_per_batch'] = max_reads_per_batch

    params['MAX_DP_BY_FIG'] = 30 # max dp for which figures plotted by each dp
    params['MAX_DP_CG_MOTIF'] = 50

    params['MAX_COORDINATE'] = 2**63-1 # default end of chrs
    params['work_dir'] = '.'
    # params['report_dir'] = os.path.join(params['work_dir'], 'report')
    params['report_dir'] = params['work_dir']
    # params['img_dir'] = os.path.join(params['report_dir'], 'img2')

    #### copy/substitue args
    for key, value in options.__dict__.items():
        params[key] = value
        if key.startswith('include_'):
            data[key] = int(value)

    # MUST FOLLOW COPYING VALUES FROM options
    data['include_mt'] = int(options.chr_MT != '-' and options.include_mt) 
    data['include_plastid'] = int(options.chr_plastid != '-' and options.include_plastid)
    data['include_lambda'] = int(options.chr_lambda != '-' and options.include_lambda)
    
    data['MAX_DP_BY_FIG'] = params['MAX_DP_BY_FIG']

    params['img_dir'] = os.path.join(params['report_dir'], 'img')
    myMakeDirs(params['work_dir'])
    myMakeDirs(params['report_dir'])
    myMakeDirs(params['img_dir'])

    data['datetime'] = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return None


class MyArgumentParser(ArgumentParser):
    def __init__(self, **args) -> None:
        super().__init__(**args)

        self.add_argument('-b', '--bam-file', dest='bamfile', help='a .bam file', type=str, required=True)
        self.add_argument('-f', '--fa-file', dest='fafile', help='a .fa[.gz] file of referece genome', type=str, required=True)
        self.add_argument('-g', '--gtf-file', dest='gtffile', help='a .gtf[.gz] file', type=str, required=False, default='-')
        self.add_argument('--chr', dest='testchrs', help='nuclear chromosomes to diagnose, seprated by comma(,), such as "chr1,chr2,chr3", all chromosomes by defaults', type=str, required=False, default='all')
        self.add_argument('--mt', dest='chr_MT', help='name of mitochondrial DNA', type=str, required=False, default='-')
        self.add_argument('--plastid', dest='chr_plastid', help='name of plastid DNA', type=str, required=False, default='-')
        self.add_argument('--control-DNA', dest='chr_lambda', help='name of spiked-in DNA, usually lambda DNA', type=str, required=False, default='-')

        self.add_argument('--diag-quality', dest='include_quality', help='diagnose sequencing and mapping quality or not, true/false, or yes/no', type=as_bool, required=False, default='yes')
        self.add_argument('--diag-pangene', dest='include_pangene', help='diagnose pangene methylation or not, true/false, or yes/no', type=as_bool, required=False, default='yes')
        self.add_argument('--diag-motif', dest='include_motif', help='diagnose CpG-motif-related patterns or not, true/false, or yes/no', type=as_bool, required=False, default='yes')
        self.add_argument('--diag-saturation', dest='include_saturation', help='diagnose sequencing saturation or not, true/false, or yes/no', type=as_bool, required=False, default='yes')
        self.add_argument('--diag-mt', dest='include_mt', help='diagnose reads mapped to mitochondrial DNA or not, true/false, or yes/no', type=as_bool, required=False, default='yes')
        self.add_argument('--diag-plastid', dest='include_plastid', help='diagnose reads mapped to plastid DNA or not, true/false, or yes/no', type=as_bool, required=False, default='yes')
        self.add_argument('--diag-control', dest='include_lambda', help='diagnose reads mapped to spiked-in control DNA (usually lambda DNA) or not, true/false, or yes/no', type=as_bool, required=False, default='yes')

        self.add_argument('--sampling-step', dest='nuclear_sampling_step', help='sampling step size of nuclear chromosomes, 1Kbp by defaults', type=int, required=False, default=1000)
        self.add_argument('--sampling-spacing', dest='nuclear_sampling_spacing', help='sampling spacing size of nuclear chromosomes, 10Kbp by defaults', type=int, required=False, default=10_000)
        self.add_argument('--bin-size', dest='binSize', help='bin size of nuclear chromosomes, 100kbp by defaults', type=int, default=100_000)
        self.add_argument('--bins-control', dest='bins_lambda', help='bins of spkied-in control DNA, 1000 bins by defaults', type=int, default=1000)
        self.add_argument('--bins-mt', dest='bins_MT', help='bins of mitochondrial DNA, 1000 bins by defaults', type=int, default=1000)
        self.add_argument('--bins-plastid', dest='bins_plastid', help='bins of plastid DNA, 1000 bins by defaults', type=int, default=1000)
        self.add_argument('--base-quality', dest='quality_threshold', help='base quality threshold', type=int, required=False, default=0)
        self.add_argument('--read-quality', dest='read_quality', help='read mapping quality threshold', type=int, required=False, default=0)
        self.add_argument('--max-depth', dest='MAXDEPTH', help='depth larger than max depth will be truncated, 200 by defaults', type=int, default=200)
        self.add_argument('--reads-for-quality', dest='reads_to_sample', help='reads to sample for quality statstics, 10K reads by defaults', type=int, default=10_000)
        self.add_argument('--num-pangene', dest='PANGENE_SAMPLED', help='genes to sample for pangene methylaion level, 1000 by defaults', type=int, default=1_000)
        self.add_argument('--max-depth-plot', dest='MAX_DP_BY_FIG', help='max depth for figures plotted for sites of each specific depth, 30 by defaults', type=int, default=20)
        self.add_argument('--max-depth-motig', dest='MAX_DP_CG_MOTIF', help='max depth in CpG-motif diagnosis, 50 by defaults', type=int, default=50)
        self.add_argument('--swap-strand', dest='swap_strand', help='swap read counts on two strands, true/false, or yes/no', type=as_bool, required=False, default='no')
        self.add_argument('-o', '--report-dir', dest='report_dir', help='report directory, "bsDoctor-report" by defaults', type=str, default='bsDoctor-report')
        # self.add_argument('--figure-subdir', dest='img_dir', help='figure subdirectory', type=str, default='img')
        self.add_argument('--ploidy', dest='ploidy', help='ploidy of the genome, 2 (diploidy) by defaults', type=int, default=2)
        self.add_argument('-v', '--version', action='version', version='%(prog)s ' + __version__)
        
        self.add_argument('--save-svg', dest='save_svg', help='save .svg figures or not, yes by defaults', type=as_bool, default='yes')

        return None

def check_args(options) -> None:
    assert not options.include_pangene or options.gtffile != '-', 'Must specify "-g/--gtf-file" to diagnose pangene methylation.'

    return None
    