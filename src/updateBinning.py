
import math
from typing import Dict

import numpy as np

from src.coverage import IntervalCoverage, CovContig, KmerCov, BinCov
from src.config import params
from src.utils import *

def update_binning_nuclear(intv: IntervalCoverage) -> None:
    dict_binning = params['dict_binning']
    MAXDEPTH = params['MAXDEPTH'] - 1
    MAX_DP_BY_FIG = params['MAX_DP_BY_FIG']

    # init for stranded CGs 
    strand_CGcov = [0,0]
    # CGid = -1
    CGmethC = -1

    key = (intv.chr, intv.bin)
    value = dict_binning.get(key)
    if value is None:
        dict_binning[key] = BinCov()
        value = dict_binning[key]
    
    value.length += intv.length
    value.dpW += intv.covWsum.sum()
    value.dpC += intv.covCsum.sum()
    depth = intv.covWsum + intv.covCsum

    # cov from 0, 
    # record dp=k instead of dp >= k
    covw = np.fmin(MAXDEPTH, intv.covWsum)
    covc = np.fmin(MAXDEPTH, intv.covCsum)
    covd = np.fmin(MAXDEPTH, intv.covWsum + intv.covCsum)
    dpmes = np.fmin(MAXDEPTH, intv.covMeth)
    # for CG strandness cov
    dpmes_strand = np.fmin(99, intv.covMeth) 
    dpmefig = np.fmin(MAX_DP_BY_FIG, intv.covMeth)

    for i in range(intv.length):
        base = intv.base[i]
        # if intv.base[i] == 'N': continue
        if base not in BASES: continue
        
        if base == 'A' or base == 'T':
            value.nAandT += 1
            value.dpAandTW += intv.covWsum[i]
            value.dpAandTC += intv.covCsum[i]
        else:
            value.nCandG += 1
            value.dpCandGW += intv.covWsum[i]
            value.dpCandGC += intv.covCsum[i]

        # value.covW[:min(MAXDEPTH, intv.covWsum[i])] += 1
        # value.covC[:min(MAXDEPTH, intv.covCsum[i])] += 1
        # value.cov[:min(MAXDEPTH, intv.covWsum[i]+intv.covCsum[i])] += 1

        # ATCG cov
        value.covW[covw[i]] += 1
        value.covC[covc[i]] += 1
        value.cov[covd[i]] += 1

        # dp = min(MAXDEPTH, intv.covMeth[i])
        # DP_FIG = min(MAX_DP_BY_FIG, intv.covMeth[i])

        dpme = dpmes[i]
        dple0 = dpme > 0
        if intv.cgContext[i] == 'CG':
            # value.nCG +=1
            # value.dpCG += intv.covMeth[i]
            # value.covnCG[dpme] += dple0
            # value.meCG[dpme] += intv.meth[i]

            if base == 'C': 
                value.nCGW += 1
                value.dpCGW += intv.covMeth[i]
                value.covnCGW[dpme] += 1
                value.meCGW[dpme] += intv.meth[i]
                
                # CG strandness
                strand_CGcov[0] = dpmes_strand[i]
                CGmethC = intv.meth[i]
            else:
                value.nCGC += 1
                value.dpCGC += intv.covMeth[i]
                value.covnCGC[dpme] += 1
                value.meCGC[dpme] += intv.meth[i]
                
                # stranded CG coverage
                strand_CGcov[1] = dpmes_strand[i]
                if CGmethC >= 0: # remove single G
                    value.stranded_CG_depth[*strand_CGcov] += 1

                # stranded CG meth
                if dple0 and min(strand_CGcov) >= 1:
                    binC = min(19, math.ceil(CGmethC / 500))
                    binG = min(19, math.ceil(intv.meth[i] / 500))
                    value.stranded_CG_meth[binC, binG] += 1
                    # CG meth diff by depth
                    md = CGmethC - intv.meth[i]
                    bindiff = min(39, math.ceil((md+10000) / 500)) # 2*20
                    value.CGmeth_diff_by_depth[bindiff, min(strand_CGcov)] += 1
                    
                strand_CGcov = [0,0] # reset
                CGmethC = -1

            # dp distrbution of low/high meth CGs
            # if not math.isnan(intv.meth[i]):
            if dple0:
                if intv.meth[i] <= 3000:
                    value.nCGLowMeth[dpme] += 1
                elif intv.meth[i] >= 7000:
                    value.nCGHighMeth[dpme] += 1
                # nCGMethBin, meth dist
                kbin = min(19, math.ceil(intv.meth[i] / 500))
                # C count by meth level bin
                value.nCGMethBin[kbin, dpmefig[i]] += 1

        elif intv.cgContext[i] == 'CHG':
            # value.nCHG +=1
            # value.dpCHG += intv.covMeth[i]
            # value.covnCHG[:dp] += 1
            # value.meCHG[:dp] += intv.meth[i]

            if base == 'C': 
                value.nCHGW += 1
                value.dpCHGW += intv.covMeth[i]
                value.covnCHGW[dpme] += 1
                value.meCHGW[dpme] += intv.meth[i]
            else:
                value.nCHGC += 1
                value.dpCHGC += intv.covMeth[i]
                value.covnCHGC[dpme] += 1
                value.meCHGC[dpme] += intv.meth[i]
            # nCGMethBin
            # if not math.isnan(intv.meth[i]):
            if dple0:
                kbin = min(19, math.ceil(intv.meth[i] / 500))
                value.nCHGMethBin[kbin, dpmefig[i]] += 1

        elif intv.cgContext[i] == 'CHH':
            # value.nCHH +=1
            # value.dpCHH += intv.covMeth[i]
            # value.covnCHH[:dp] += 1
            # value.meCHH[:dp] += intv.meth[i]

            if base == 'C': 
                value.nCHHW += 1
                value.dpCHHW += intv.covMeth[i]
                value.covnCHHW[dpme] += 1
                value.meCHHW[dpme] += intv.meth[i]
            else:
                value.nCHHC += 1
                value.dpCHHC += intv.covMeth[i]
                value.covnCHHC[dpme] += 1
                value.meCHHC[dpme] += intv.meth[i]
            # nCGMethBin
            # if not math.isnan(intv.meth[i]):
            if dple0:
                kbin = min(19, math.ceil(intv.meth[i] / 500))
                value.nCHHMethBin[kbin, dpmefig[i]] += 1

        # error bases for A/T sites (DP>=5)
        elif depth[i] >= 5:
            value.ATdp += depth[i]
            # value.misbase += depth - max(max(intv.covW[:,i]), max(intv.covC[:,i]))
            value.misbase += depth[i] - np.max(intv.covW[:,i] + intv.covC[:,i])
    return None


# def update_binning_contig(intv: IntervalCoverage) -> None:
#     dict_binning = params['dict_binning']
#     MAXDEPTH = params['MAXDEPTH']
#     MAX_DP_BY_FIG = params['MAX_DP_BY_FIG'] + 1

#     # # init for stranded CGs 
#     # strand_CGcov = [0,0]
#     # CGid = -1
#     # CGmethC = -1

#     for i in range(intv.length):
#         if intv.base[i] == 'N': continue
#         # if intv.covMeth[i] == 0: continue

#         key = (intv.chr, intv.bin[i])
#         value = dict_binning.get(key)
#         if value is None:
#             dict_binning[key] = BinCov()
#             value = dict_binning[key]
            
#         value.length += 1
#         depth = intv.covWsum[i] + intv.covCsum[i]
#         value.dp += depth
#         value.dpW += intv.covWsum[i]
#         value.dpC += intv.covCsum[i]

#         value.covW[:min(MAXDEPTH, intv.covWsum[i])] += 1
#         value.covC[:min(MAXDEPTH, intv.covCsum[i])] += 1
#         value.cov[:min(MAXDEPTH, intv.covWsum[i]+intv.covCsum[i])] += 1
        
#         dp = min(MAXDEPTH, intv.covMeth[i])
#         DP_FIG = min(MAX_DP_BY_FIG, intv.covMeth[i])

#         if intv.cgContext[i] == 'CG':
#             value.nCG +=1
#             value.dpCG += intv.covMeth[i]
#             value.covnCG[:dp] += 1
#             value.meCG[:dp] += intv.meth[i] # if dp==0, nothing to do

#             if intv.base[i] == 'C': 
#                 value.nCGW += 1
#                 value.dpCGW += intv.covMeth[i]
#                 value.covnCGW[:dp] += 1
#                 value.meCGW[:dp] += intv.meth[i]
                
#                 # # CG strandness
#                 # strand_CGcov[0] = min(99, intv.covMeth[i])
#                 # CGid = intv.start + i
#                 # CGmethC = intv.meth[i]
#             else:
#                 value.nCGC += 1
#                 value.dpCGC += intv.covMeth[i]
#                 value.covnCGC[:dp] += 1
#                 value.meCGC[:dp] += intv.meth[i]
                
#                 # # stranded CG coverage
#                 # strand_CGcov[1] = min(99, intv.covMeth[i])
#                 # value.stranded_CG_depth[*strand_CGcov] += 1

#                 # # stranded CG meth
#                 # if min(strand_CGcov) >= 1:
#                 #     binC = min(19, math.ceil(CGmethC * 20))
#                 #     binG = min(19, math.ceil(intv.meth[i] * 20))
#                 #     value.stranded_CG_meth[binC, binG] += 1
#                 #     # CG meth diff by depth
#                 #     md = CGmethC - intv.meth[i]
#                 #     bindiff = min(39, math.ceil((md+1) * 20)) # 2*20
#                 #     value.CGmeth_diff_by_depth[bindiff, min(strand_CGcov)] += 1
                                    
#                 # strand_CGcov = [0,0] # reset

#             # # dp distrbution of low/high meth CGs
#             # if not math.isnan(intv.meth[i]):
#             #     if intv.meth[i] <= 0.3:
#             #         value.nCGLowMeth[dp - 1] += 1
#             #     elif intv.meth[i] >= 0.7:
#             #         value.nCGHighMeth[dp - 1] += 1
#             #     # nCGMethBin
#             #     kbin = min(19, math.ceil(intv.meth[i] * 20))
#             #     value.nCGMethBin[kbin, DP_FIG-1] += 1

#         elif intv.cgContext[i] == 'CHG':
#             value.nCHG +=1
#             value.dpCHG += intv.covMeth[i]
#             value.covnCHG[:dp] += 1
#             value.meCHG[:dp] += intv.meth[i]

#             if intv.base[i] == 'C': 
#                 value.nCHGW += 1
#                 value.dpCHGW += intv.covMeth[i]
#                 value.covnCHGW[:dp] += 1
#                 value.meCHGW[:dp] += intv.meth[i]
#             else:
#                 value.nCHGC += 1
#                 value.dpCHGC += intv.covMeth[i]
#                 value.covnCHGC[:dp] += 1
#                 value.meCHGC[:dp] += intv.meth[i]
#             # nCGMethBin
#             if not math.isnan(intv.meth[i]):
#                 kbin = min(19, math.ceil(intv.meth[i] * 20))
#                 value.nCHGMethBin[kbin, DP_FIG-1] += 1
#         elif intv.cgContext[i] == 'CHH':
#             value.nCHH +=1
#             value.dpCHH += intv.covMeth[i]
#             value.covnCHH[:dp] += 1
#             value.meCHH[:dp] += intv.meth[i]

#             if intv.base[i] == 'C': 
#                 value.nCHHW += 1
#                 value.dpCHHW += intv.covMeth[i]
#                 value.covnCHHW[:dp] += 1
#                 value.meCHHW[:dp] += intv.meth[i]
#             else:
#                 value.nCHHC += 1
#                 value.dpCHHC += intv.covMeth[i]
#                 value.covnCHHC[:dp] += 1
#                 value.meCHHC[:dp] += intv.meth[i]
#             # nCGMethBin
#             if not math.isnan(intv.meth[i]):
#                 kbin = min(19, math.ceil(intv.meth[i] * 20))
#                 value.nCHHMethBin[kbin, DP_FIG-1] += 1

#         # error bases for A/T sites (DP>=5)
#         elif depth >= 5:
#             value.ATdp += depth
#             # value.misbase += depth - max(max(intv.covW[:,i]), max(intv.covC[:,i]))
#             value.misbase += depth - max(intv.covW[:,i] + intv.covC[:,i])
#     return None

def update_cgkmer(intv: IntervalCoverage) -> None:
    dict_cgkmer = params['dict_cgkmer']
    MAXDEPTH = params['MAXDEPTH'] - 1
    
    dps = np.fmin(MAXDEPTH, intv.covMeth)

    for i in range(intv.length):
        # if intv.cgContext[i] != 'CG': continue
        if intv.cgkmer[i] == '-': continue
        dp = dps[i]
        dple0 = dp > 0

        key = intv.cgkmer[i]
        value = dict_cgkmer.get(key)
        if value is None:
            dict_cgkmer[key] = KmerCov()
            value = dict_cgkmer[key]
            
        # DP = np.arange(start=1, stop=101)
        # MAXDP = 100
        # dp = min(MAXDEPTH, intv.covMeth[i])
        # value.n += 1
        # value.dp += intv.covMeth[i]
        # value.cov[:dp] += 1
        # value.me[:dp] += intv.meth[i]
        if intv.base[i] == 'C':
            value.nW += 1
            value.dpW += intv.covMeth[i]
            value.covW[dp] += 1
            value.meW[dp] += intv.meth[i]
        else:
            value.nC += 1
            value.dpC += intv.covMeth[i]
            value.covC[dp] += 1
            value.meC[dp] += intv.meth[i]
    return None

def update_contig(cov: CovContig, intv: IntervalCoverage) -> None:
    for i in range(intv.length):
        # if intv.base[i] == 'N': continue

        # bin index
        kbin = intv.bin[i]

        cov.length[kbin] += 1
        dpWC = intv.covWsum[i] + intv.covCsum[i]
        cov.dp[kbin] += dpWC
        cov.dpW[kbin] += intv.covWsum[i]
        cov.dpC[kbin] += intv.covCsum[i]
        if dpWC > 0: cov.cov[kbin] += 1

        # base-wise error rate
        if dpWC >= 20:
            cov.dp20[kbin] += dpWC
            cov.misbase[kbin] += dpWC - np.max(intv.covW[:,i] + intv.covC[:,i])

        if intv.base[i] == 'C' or intv.base[i] == 'G':
            cov.nC +=1
            if intv.covMeth[i] >= 1:
                cov.covnC += 1
                cov.meC += intv.meth[i]
    return None
