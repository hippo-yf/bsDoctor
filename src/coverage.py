
import math
import gzip
import re

from typing import NamedTuple, List, Tuple, Dict, Any
import pysam
import numpy as np

from src.utils import *
from src.config import params
# from src.quality import Quality

################################################
## Quality
################################################

class Quality():
    def __init__(self) -> None:
        self.base_quality = list()
        self.read_length = list()
        self.nreads: int = 0
        self.ndup: int = 0
        self.nqcfail: int = 0
        self.npaired: int = 0
        self.nproperpaired: int = 0
        # self.nbases: int = 0
        self.ncigars: NDArray = np.zeros((10,), dtype=int64)
        self.map_quality = list()
        return None
    
    def regularization(self) -> None:
        self.base_quality = np.hstack(self.base_quality, dtype=int16)
        self.map_quality = np.hstack(self.map_quality, dtype=int16)
        self.read_length = np.array(self.read_length, dtype=int32)
        return None


################################################
## fa methods
################################################

# interval must be involved in single chr
class GenomicInterval(NamedTuple):
    chr: str
    chr_length: int
    start: int
    end: int

# `ref` includes ref bases from `start`-2 to `end`-2
class FaGenomicInterval(NamedTuple):
    chr: str
    start: int
    end: int
    bases: str

class MyFastaFile(pysam.FastaFile):

    def rich_fetch(self, intrv: GenomicInterval, padding:int) -> str:

        bases = self.fetch(reference=intrv.chr, 
                           start=max(0, intrv.start-padding), 
                           end=min(intrv.chr_length, intrv.end+padding)
                           ).upper()
        
        # padding N
        if intrv.start < padding:
            bases = "N"*(padding-intrv.start) + bases
        if intrv.end+padding > intrv.chr_length:
            bases = bases + "N"*(intrv.end+padding-intrv.chr_length)

        return bases

################################################
## genome sampling
################################################

class GenomicIntervalGenerator:
    
    def __init__(self, 
                 fa: pysam.FastaFile, 
                 chrs,
                 start: int,
                 end: int,
                 step: int,
                 spacing: int = 0
                 ) -> None:
        
        self.chrs = fa.references
        self.lens = fa.lengths
        # if chrs == "all":
        #     self.chrs_selected = list(self.chrs)
        # else:
        #     self.chrs_selected = chrs.split(',')
        # self.chrs_selected = params['chrs_valid']
        self.chrs_selected = chrs

        self.start = start
        self.end = end
        self.step = step
        self.spacing = spacing

        assert(step > 0 and start < end)
        assert(len(self.chrs) > 0 and len(self.chrs) == len(self.lens))
        assert(isinstance(spacing, int) and spacing >= 0)
        return None
    
    def __iter__(self):
        for (chr, len) in zip(self.chrs, self.lens):
            if chr not in self.chrs_selected:
                continue

            end2 = min(self.end, len)
            start = self.start + GEN.integers(low=0, high=self.spacing*2+1, size=1)[0]
            end = start + self.step
            while start < end2:
                # if start < end2:
                end = min(end, end2)
                yield GenomicInterval(chr=chr, chr_length=len, start=start, end=end)
                start = end + GEN.integers(low=0, high=self.spacing*2+1, size=1)[0]
                end = start + self.step
                # else:
                #     break

# sample batches within bin
def cut(start: int, end: int, binSize: int):
    assert(start < end and binSize > 0)
    yield start
    i = math.ceil((start+1) / binSize)
    while(i*binSize <= end):
        if 0 < i*binSize - start < binSize/2:
            i += 1
        elif end < i*binSize + binSize/2:
            i += 1
        else:
            yield i*binSize
            i += 1
    yield end


class GenomicIntervalGeneratorWithinBin:
    def __init__(self, 
                 fa: pysam.FastaFile, 
                 chrs,
                 start: int,
                 end: int,
                 step: int,
                 spacing: int,
                 binSize: int
                 ) -> None:
        
        self.chrs = fa.references
        self.lens = fa.lengths
        self.chrs_selected = chrs

        self.start = start
        self.end = end
        self.step = step
        self.spacing = spacing
        self.binSize = binSize
        assert(step > 0 and start < end and binSize > 0)
        assert(len(self.chrs) > 0 and len(self.chrs) == len(self.lens))
        assert(isinstance(spacing, int) and spacing >= 0)
        return None
    
    def __iter__(self):
        for (contig, llen) in zip(self.chrs, self.lens):
            if contig not in self.chrs_selected:
                continue
            cuts = list(cut(self.start, min(self.end, llen), self.binSize))
            for ib in range(len(cuts) - 1):
                end2 = cuts[ib+1]
                start = cuts[ib] + GEN.integers(low=0, high=self.spacing*2+1, size=1)[0]
                end = start + self.step
                while start < end2:
                    end = min(end, end2)
                    yield GenomicInterval(chr=contig, chr_length=llen, start=start, end=end)
                    start = end + GEN.integers(low=0, high=self.spacing*2+1, size=1)[0]
                    end = start + self.step



## sampling for saturation curve 
def randomGenomicRegion(fa: pysam.FastaFile, bam: pysam.AlignmentFile, size=100):
    chrs = fa.references
    lens = fa.lengths
    reads = np.zeros((len(chrs),), dtype=int)
    mappedreads = {}
    for c in bam.get_index_statistics():
        mappedreads[c.contig] = c.mapped
    for i, c in enumerate(chrs):
        value = mappedreads.get(c)
        if value is None: reads[i] = 0
        else: reads[i] = value
    # prop
    p = reads/np.sum(reads)
    idx = list(range(len(chrs)))
    while(True):
        ci = random.choices(idx, weights=p, k=1)[0]
        if lens[ci] <= size:
            a, b = 0, lens[ci]
        else:
            pos = random.randrange(0, lens[ci])
            a, b = max(0, pos-int(size/2)), min(lens[ci], pos+int(size/2))
        yield GenomicInterval(chr=chrs[ci], chr_length=lens[ci], start=a, end=b)
    
################################################
## gene samlping for pangene meth
################################################

## gene intervals
class Gene():
    def __init__(self, chr, chr_length, start, end, strand, biotype=None) -> None:
        self.chr: str = chr
        self.chr_length: int = chr_length
        self.start: int = start
        self.end: int = end
        self.strand: str = strand
        self.biotype: str | None = biotype
        
        # 
        gene_padding = params['gene_padding']
        gene_breaks = params['gene_breaks']

        self.start_padding = max(0, start-gene_padding)
        self.end_padding = min(chr_length, end+gene_padding)
        breaks1 = np.linspace(self.start_padding-1, start, num=gene_breaks+1)
        breaks2 = np.linspace(start, end, num=gene_breaks+1)
        breaks3 = np.linspace(end, self.end_padding+1, num=gene_breaks+1)
        self.breaks = np.block([breaks1, breaks2[1:], breaks3[1:]])
        # bins = np.zeros(gene_breaks*3, dtype=int)
        if strand == '+':
            bins = np.arange(0, gene_breaks*3, dtype=int)
        else:
            bins = -np.arange(-(gene_breaks*3-1), 1)
        self.bins = bins
        return None


# parse gtf file
def genesGenerator(file: str):
    fa = params['fa']

    gtf = gzip.open(file, 'rt')
    while line := gtf.readline().strip():
        if line.startswith('#'): continue
        keys = line.split('\t')
        if keys[2] != 'gene': continue
        # m = re.search(r'gene_type=([^;]*);', keys[8]) # gff3
        m = re.search(r'gene_(?:bio)?type "([^;]*)";', keys[8]) # gtf
        if m is None: biotype = None
        else: biotype = m.group(1)
        chr_length = fa.get_reference_length(keys[0])
        yield Gene(chr=keys[0], chr_length=chr_length, start=int(keys[3]), end=int(keys[4]), strand=keys[6], biotype=biotype)
    gtf.close()


## coverage of pangene
class CovPanGene():
    def __init__(self, bins: int) -> None:
        # bins = params['gene_breaks']*3
        
        self.nCG: NDArray[int32] = np.zeros((bins,), dtype=int32)
        self.meCG: NDArray[int64] = np.zeros((bins,), dtype=int64)
        self.nCGW: NDArray[int32] = np.zeros((bins,), dtype=int32)
        self.meCGW: NDArray[int64] = np.zeros((bins,), dtype=int64)
        self.nCGC: NDArray[int32] = np.zeros((bins,), dtype=int32)
        self.meCGC: NDArray[int64] = np.zeros((bins,), dtype=int64)
        self.nCHG: NDArray[int32] = np.zeros((bins,), dtype=int32)
        self.meCHG: NDArray[int64] = np.zeros((bins,), dtype=int64)
        self.nCHGW: NDArray[int32] = np.zeros((bins,), dtype=int32)
        self.meCHGW: NDArray[int64] = np.zeros((bins,), dtype=int64)
        self.nCHGC: NDArray[int32] = np.zeros((bins,), dtype=int32)
        self.meCHGC: NDArray[int64] = np.zeros((bins,), dtype=int64)
        self.nCHH: NDArray[int32] = np.zeros((bins,), dtype=int32)
        self.meCHH: NDArray[int64] = np.zeros((bins,), dtype=int64)
        self.nCHHW: NDArray[int32] = np.zeros((bins,), dtype=int32)
        self.meCHHW: NDArray[int64] = np.zeros((bins,), dtype=int64)
        self.nCHHC: NDArray[int32] = np.zeros((bins,), dtype=int32)
        self.meCHHC: NDArray[int64] = np.zeros((bins,), dtype=int64)
        return None
        
# cov_pangene = CovPanGene()

class KmerCov():
    def __init__(self) -> None:
        MAXDEPTH = params['MAXDEPTH']

        # self.n: int = 0
        self.nW: int = 0
        self.nC: int = 0
        # self.dp: int = 0
        self.dpW: int = 0
        self.dpC: int = 0

        # init = np.array([0]*MAXDEPTH, dtype=np.int64)
        # init2 = np.array([0]*MAXDEPTH, dtype=np.float64)
        # self.cov: NDArray[np.int64] = init.copy()
        self.covW: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        self.covC: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        # self.me: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        self.meW: NDArray[int64] = np.zeros((MAXDEPTH,), dtype=int64)
        self.meC: NDArray[int64] = np.zeros((MAXDEPTH,), dtype=int64)
        return None
    
class BinCov():
    def __init__(self) -> None:  
        MAXDEPTH = params['MAXDEPTH']
        MAX_DP_BY_FIG = params['MAX_DP_BY_FIG'] + 1

        # init = np.zeros((MAXDEPTH,), dtype=np.int64)
        # init2 = np.zeros((MAXDEPTH,), dtype=np.int16)
 
        self.length: int = 0 # bases sampled
        # self.nCG: int = 0 # CGs/CHGs/CHHs in the genome
        self.nCGW: int = 0
        self.nCGC: int = 0
        # self.nCHG: int = 0 # CGs/CHGs/CHHs in the genome
        self.nCHGW: int = 0
        self.nCHGC: int = 0
        # self.nCHH: int = 0 # CGs/CHGs/CHHs in the genome
        self.nCHHW: int = 0
        self.nCHHC: int = 0

        # self.dp: int = 0 # total DP
        self.dpW: int = 0
        self.dpC: int = 0
        # self.dpCG: int = 0 # totol CG dp
        self.dpCGW: int = 0
        self.dpCGC: int = 0
        # self.dpCHG: int = 0 # totol CHG dp
        self.dpCHGW: int = 0
        self.dpCHGC: int = 0
        # self.dpCHH: int = 0 # totol CHG dp
        self.dpCHHW: int = 0
        self.dpCHHC: int = 0
        
        # ATCG coverage
        self.cov: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        self.covW: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        self.covC: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)

        # CGs coverage and meth
        # self.covnCG: NDArray[np.int64] = init.copy()
        self.covnCGW: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        self.covnCGC: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        # self.meCG: NDArray[np.int64] = init.copy()
        self.meCGW: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        self.meCGC: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)

        # CHGs coverage and meth
        # self.covnCHG: NDArray[np.int64] = init.copy()
        self.covnCHGW: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        self.covnCHGC: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        # self.meCHG: NDArray[np.int64] = init.copy()
        self.meCHGW: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        self.meCHGC: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)

        # CHHs coverage and meth
        # self.covnCHH: NDArray[np.int64] = init.copy()
        self.covnCHHW: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        self.covnCHHC: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        # self.meCHH: NDArray[np.int64] = init.copy()
        self.meCHHW: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        self.meCHHC: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)

        # coverage of low/high meth CGs
        self.nCGLowMeth: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)
        self.nCGHighMeth: NDArray[int32] = np.zeros((MAXDEPTH,), dtype=int32)

        # strand-specific coverage and methylation
        self.stranded_CG_depth: NDArray[int32] = np.zeros((100,100), dtype=int32)
        # 20X20 bins of [0,1]X[0,1]
        self.stranded_CG_meth: NDArray[int32] = np.zeros((20,20), dtype=int32)
        # 40(+-) row meth diff bins X 100 cols of DP
        self.CGmeth_diff_by_depth: NDArray[int32] = np.zeros((40,100), dtype=int32)

        # misbases for ref A/T
        self.misbase: int = 0 # inconsistent bases
        self.ATdp: int = 0 # total depths/bases of A/T

        # meth distribution, 20 meth bins of [0,1] X depth [1,30]
        self.nCGMethBin: NDArray[int32] = np.zeros((20, MAX_DP_BY_FIG), dtype=int32)
        self.nCHGMethBin: NDArray[int32]  = np.zeros((20, MAX_DP_BY_FIG), dtype=int32)
        self.nCHHMethBin: NDArray[int32]  = np.zeros((20, MAX_DP_BY_FIG), dtype=int32)
        return None

# for whole contig not bin
class CovContig():
    def __init__(self, contig:str) -> None:
        binsContig = params['binsContig']

        bins = binsContig[contig]
        # len of each bin
        self.length: NDArray[int32] = np.zeros((bins,), dtype=int32) 
        # covered
        self.cov: NDArray[int32] = np.zeros((bins,), dtype=int32)
        # Cytosines in total
        self.nC: int = 0 

        # total dp for each bin
        self.dp: NDArray[int64] = np.zeros((bins,), dtype=int64)
        self.dpW: NDArray[int64] = np.zeros((bins,), dtype=int64)
        self.dpC: NDArray[int64] = np.zeros((bins,), dtype=int64)

        # Cs coverage and meth
        self.covnC: int = 0
        self.meC: int64 = int64(0)

        ## mistake bases, different with ref/max base
        # nbases with DP>=20
        self.dp20: NDArray[int32] = np.zeros((bins,), dtype=int32) 
        self.misbase: NDArray[int64] = np.zeros((bins,), dtype=int64)
        return None

##############################################
## bam methods
##############################################

class Coverage(NamedTuple):
    watson: NDArray[int32]
    crick: NDArray[int32]

def check_read(forward_read: bool, read_quality: int):
    # def valid_read(read: pysam.AlignedSegment):
    #     return (forward_read ^ read.is_reverse) and (not read.is_unmapped) and (not read.is_duplicate) and (not read.is_secondary) and (not read.is_qcfail) and (read.mapping_quality>=read_quality)
    def valid_read(read: pysam.AlignedSegment):
        return (forward_read ^ read.is_reverse) and (read.mapping_quality>=read_quality)
    return valid_read

class IntervalCoverage(NamedTuple):
    chr: str
    length: int
    start: int
    end: int
    bin: NDArray[int32] | int
    base: str
    cgContext: List[str]
    cgkmer: List[str]
    covW: NDArray[int32]
    covC: NDArray[int32]
    covWsum: NDArray[int32]
    covCsum: NDArray[int32]
    covMeth: NDArray[int32]
    meth: NDArray[int16]

    def printInterval(self) -> None:
        for i in range(self.length):
            if len(bin) == 1:
                ib = bin
            else:
                ib = bin[i]
            print(f'{self.chr}\t{self.start+i}\t{ib}\t{self.base[i]}\t{self.cgContext[i]}\t{self.cgkmer[i]}\t{self.covWsum[i]}\t{self.covCsum[i]}\t{self.covMeth[i]}\t{self.meth[i]:.2f}')
        return None


class MyAlignmentFile(pysam.AlignmentFile):

    def Watson_Crick_coverage(self, intrv: GenomicInterval) -> Coverage:
        quality_threshold = params['quality_threshold']
        read_quality = params['read_quality']
        swap_strand = params['swap_strand']

        # 4 X nbases
        cov_watson = self.count_coverage(
            contig=intrv.chr, 
            start=intrv.start, 
            stop=intrv.end,
            quality_threshold=quality_threshold,
            read_callback=check_read(forward_read=True, read_quality=read_quality)
            )
        cov_crick = self.count_coverage(
            contig=intrv.chr, 
            start=intrv.start, 
            stop=intrv.end, 
            quality_threshold=quality_threshold,
            read_callback=check_read(forward_read=False, read_quality=read_quality)
            )
        
        # in some bams, the read strandness seem be reversly flaged in `FLAG`
        # parsed in read.is_reverse

        if swap_strand:
            cov_watson, cov_crick = cov_crick, cov_watson
        return Coverage(
            watson=np.array(cov_watson, dtype=np.int32),
            crick=np.array(cov_crick, dtype=np.int32)
            )
    
    def detailedCoverage(self, interval: GenomicInterval) -> IntervalCoverage:
        fa = params['fa']
        binSize = params['binSize']
        cgkmerSize = params['cgkmerSize']

        # context size
        consize = params['context_size']
        # ref sequences
        bases = fa.rich_fetch(interval, padding=cgkmerSize//2).upper()

        # bam coverages
        covs = self.Watson_Crick_coverage(interval)
        cov_sum_W = np.sum(covs.watson, axis=0, dtype=np.int32)
        cov_sum_C = np.sum(covs.crick, axis=0, dtype=np.int32)

        CG_contexts = []
        cgkmers = []
        length = interval.end - interval.start
        covMeths = np.zeros((length,), dtype=np.int32)
        # meth = np.array([math.nan]*length, dtype=np.float32)
        meth = np.zeros((length,), dtype=np.int16)
        # bin = np.zeros((length,), dtype=np.int32)

        # bin index
        # kbin = interval.start // binSize
        nbases = 0 # bases taken
        # keep = np.array([False]*length, dtype=np.bool_) # which bases to take
        
        for i in range(length):
            # skip A/Ts
            j = i + cgkmerSize//2
            base = bases[j]
            # if base == 'A' or base == 'T' or base == 'N': continue

            # # bin index
            # if i + interval.start >= (kbin + 1) * binSize:
            #     kbin += 1
            # bin[i] = kbin

            CG_context = "-"
            # meth_ratio = math.nan
            cgkmer = '-'
            nCT = 0
            
            # if base == 'C' or base == 'G':
            if base == 'C':
                # CG/CHG/CHH
                bases_con = bases[j:(j+consize)]
                # CG_context = '-' if 'N' in bases_con else CG_CONTEXT_FORWARD_HASH[bases_con]
                CG_context = CG_CONTEXT_FORWARD_HASH[bases_con] if all([b in BASES for b in bases_con]) else '-' 

                # dinucleatide context CA/CT/...
                # dicontext = bases[j:(j+2)]
                nCT = covs.watson[1,i] + covs.watson[3,i] # order of ACGT
                nC = covs.watson[1,i]
                if nCT > 0: meth[i] = np.int16(nC*10_000/nCT )
                # CG kmer
                if CG_context == 'CG':
                    cgkmer = bases[(j-2):(j+4)]
                    if 'N' in cgkmer: cgkmer = '-'
            elif base == 'G':
                bases_con = bases[(j-consize+1):(j+1)]
                # CG_context = '-' if 'N' in bases_con else CG_CONTEXT_REVERSE_HASH[bases_con]
                CG_context = CG_CONTEXT_REVERSE_HASH[bases_con] if all([b in BASES for b in bases_con]) else '-' 
                nC = covs.crick[2,i]
                nCT = covs.crick[0,i] + covs.crick[2,i]
                if nCT > 0: meth[i] = np.int16(nC*10_000/nCT )
                if CG_context == 'CG':
                    kmer = bases[(j-3):(j+3)]
                    if 'N' in kmer: 
                        cgkmer = '-'
                    else:
                        cgkmer = reverseComp(kmer)
            cgkmers.append(cgkmer)
            CG_contexts.append(CG_context)
            covMeths[i] = nCT

        return IntervalCoverage(
            chr=interval.chr, 
            length=length,
            start=interval.start, 
            end=interval.end, 
            bin = interval.start // binSize,
            base=bases[(cgkmerSize//2):(length+cgkmerSize//2)],
            cgContext=CG_contexts, 
            cgkmer=cgkmers, 
            covW=covs.watson,
            covC=covs.crick,
            covWsum=cov_sum_W, 
            covCsum=cov_sum_C, 
            covMeth=covMeths,
            meth=meth
            )
    
    # for extranuclear contigs: MT, plastid, lambda DNA
    def detailedCoverageContig(self, interval: GenomicInterval) -> IntervalCoverage:
        fa = params['fa']
        binSizeContig = params['binSizeContig']
        cgkmerSize = params['cgkmerSize']

        # ref sequences
        bases = fa.rich_fetch(interval, padding=cgkmerSize//2).upper()

        # bam coverages
        covs = self.Watson_Crick_coverage(interval)
        cov_sum_W = np.sum(covs.watson, axis=0, dtype=int32)
        cov_sum_C = np.sum(covs.crick, axis=0, dtype=int32)

        length = interval.end - interval.start
        covMeths = np.zeros((length,), dtype=int32)
        meth = np.zeros((length,), dtype=int16)
        bin = np.zeros((length,), dtype=int32)

        # bin index 
        # binSize = binSizeContig[chrs_alias[interval.chr]]
        binSize = binSizeContig[interval.chr]

        kbin = interval.start // binSize
        
        for i in range(length):
            j = i + cgkmerSize//2
            base = bases[j]
            # if base == 'N': continue # processed in next step

            # bin index
            if i + interval.start >= (kbin + 1) * binSize:
                kbin += 1
            bin[i] = kbin

            if base == 'C' or base == 'G':
                # meth_ratio = math.nan
                nCT = 0                
                if base == 'C':
                    # dinucleatide context CA/CT/...
                    # dicontext = bases[j:(j+2)]
                    nCT = covs.watson[1,i] + covs.watson[3,i]
                    nC = covs.watson[1,i]
                else:
                    nC = covs.crick[2,i]
                    nCT = covs.crick[0,i] + covs.crick[2,i]
                meth_ratio = np.int16(nC/nCT*10_000) if nCT > 0 else 0 
                meth[i] = meth_ratio
                covMeths[i] = nCT

        return IntervalCoverage(
            chr=interval.chr, 
            length=length,
            start=interval.start, 
            end=interval.end, 
            bin = bin,
            base=bases[(cgkmerSize//2):(length+cgkmerSize//2)],
            cgContext=[''], 
            cgkmer=[''], 
            covW=covs.watson,
            covC=covs.crick,
            covWsum=cov_sum_W, 
            covCsum=cov_sum_C, 
            covMeth=covMeths,
            meth=meth
            )

    def update_pangene(self, cov: CovPanGene, g: Gene) -> None:
        # context size
        consize = params['context_size']
        fa = params['fa']
        cgkmerSize = params['cgkmerSize']

        intv = GenomicInterval(chr=g.chr, chr_length=g.chr_length, start=g.start_padding, end=g.end_padding)
        
        # ref sequences
        bases = fa.rich_fetch(intv, padding=cgkmerSize//2).upper()

        # bam coverages
        covs = self.Watson_Crick_coverage(intv)
        
        for i in range(g.end_padding - g.start_padding):
            j = i + cgkmerSize//2

            base = bases[j]
            if base == 'C' or base == 'G':
                if base == 'C':
                    # CG/CHG/CHH
                    # only CG for pangene
                    bases_con = bases[j:(j+consize)]
                    CG_context = CG_CONTEXT_FORWARD_HASH[bases_con] if all([n in BASES for n in bases_con]) else '-' 
                    # if CG_context != 'CG': continue
                    nCT = covs.watson[1,i] + covs.watson[3,i]
                    nC = covs.watson[1,i]
                else:
                    bases_con = bases[(j-consize+1):(j+1)]
                    CG_context = CG_CONTEXT_REVERSE_HASH[bases_con] if all([n in BASES for n in bases_con]) else '-' 
                    # if CG_context != 'CG': continue
                    nC = covs.crick[2,i]
                    nCT = covs.crick[0,i] + covs.crick[2,i]

                # pdb.set_trace()
                if nCT <= 0: continue
                meth_ratio = int64(nC*10000 / nCT)
                kbin = g.bins[np.argmax(g.start_padding + i < g.breaks) - 1]
                if CG_context == 'CG':
                    cov.nCG[kbin] += 1
                    cov.meCG[kbin] += meth_ratio
                    if base == 'C':
                        cov.nCGW[kbin] += 1
                        cov.meCGW[kbin] += meth_ratio
                    else:
                        cov.nCGC[kbin] += 1
                        cov.meCGC[kbin] += meth_ratio
                elif CG_context == 'CHG':
                    cov.nCHG[kbin] += 1
                    cov.meCHG[kbin] += meth_ratio
                    if base == 'C':
                        cov.nCHGW[kbin] += 1
                        cov.meCHGW[kbin] += meth_ratio
                    else:
                        cov.nCHGC[kbin] += 1
                        cov.meCHGC[kbin] += meth_ratio
                else:
                    cov.nCHH[kbin] += 1
                    cov.meCHH[kbin] += meth_ratio
                    if base == 'C':
                        cov.nCHHW[kbin] += 1
                        cov.meCHHW[kbin] += meth_ratio
                    else:
                        cov.nCHHC[kbin] += 1
                        cov.meCHHC[kbin] += meth_ratio      
        return None


    # for reads and bases sampling
    def sample_reads(self, quality: Quality, gi: GenomicInterval) -> None:
        max_reads_per_batch = params['max_reads_per_batch']

        reads = self.fetch(contig=gi.chr, start=gi.start, stop=gi.end)
        k = 0
        for r in reads:
            if k > max_reads_per_batch: return
            k += 1
            quality.nreads += 1
            quality.nqcfail += r.is_qcfail
            quality.npaired += r.is_paired
            quality.nproperpaired += r.is_proper_pair
            quality.ndup += r.is_duplicate

            quality.read_length.append(r.query_length)
            quality.base_quality.append(r.query_qualities)
            quality.map_quality.append(r.mapping_quality)

            # cigar
            cigar = r.cigartuples
            if cigar is not None:
                for c in cigar:
                    if c[0] < 10:
                        quality.ncigars[c[0]] += c[1]
        return None
