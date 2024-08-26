
import math
from typing import Dict

from src.coverage import IntervalCoverage, CovLambda, KmerCov, BinCov
from src.config import params


def update_binning(intv: IntervalCoverage):
    dict_binning = params['dict_binning']
    MAXDEPTH = params['MAXDEPTH']

    # init for stranded CGs 
    strand_CGcov = [0,0]
    CGid = -1
    CGmethC = -1

    for i in range(intv.length):
        if intv.base[i] == 'N': continue
        # if intv.covMeth[i] == 0: continue

        key = (intv.chr, intv.bin[i])
        value = dict_binning.get(key)
        if value is None:
            dict_binning[key] = BinCov()
            value = dict_binning[key]
            
        value.length += 1
        depth = intv.covWsum[i] + intv.covCsum[i]
        value.dp += depth
        value.dpW += intv.covWsum[i]
        value.dpC += intv.covCsum[i]

        # DP = np.arange(start=1, stop=101)
        # MAXDP = 100
        value.covW[:min(MAXDEPTH, intv.covWsum[i])] += 1
        value.covC[:min(MAXDEPTH, intv.covCsum[i])] += 1
        value.cov[:min(MAXDEPTH, intv.covWsum[i]+intv.covCsum[i])] += 1
        
        dp = min(MAXDEPTH, intv.covMeth[i])
        if intv.cgContext[i] == 'CG':
            value.nCG +=1
            value.dpCG += intv.covMeth[i]
            value.covnCG[:dp] += 1
            value.meCG[:dp] += intv.meth[i] # if dp==0, nothing to do

            if intv.base[i] == 'C': 
                value.nCGW += 1
                value.dpCGW += intv.covMeth[i]
                value.covnCGW[:dp] += 1
                value.meCGW[:dp] += intv.meth[i]
                
                # CG strandness
                strand_CGcov[0] = min(99, intv.covMeth[i])
                CGid = intv.start + i
                CGmethC = intv.meth[i]
            else:
                value.nCGC += 1
                value.dpCGC += intv.covMeth[i]
                value.covnCGC[:dp] += 1
                value.meCGC[:dp] += intv.meth[i]
                
                # stranded CG coverage
                strand_CGcov[1] = min(99, intv.covMeth[i])
                value.stranded_CG_depth[*strand_CGcov] += 1

                # stranded CG meth
                # if (intv.start + i == CGid + 1) and (not math.isnan(CGmethC)) and (not math.isnan(intv.meth[i])):
                if min(strand_CGcov) >= 1:
                    binC = min(19, math.ceil(CGmethC * 20))
                    binG = min(19, math.ceil(intv.meth[i] * 20))
                    value.stranded_CG_meth[binC, binG] += 1
                    # CG meth diff by depth
                    md = CGmethC - intv.meth[i]
                    bindiff = min(39, math.ceil((md+1) * 20)) # 2*20
                    value.CGmeth_diff_by_depth[bindiff, min(strand_CGcov)] += 1
                
                    
                strand_CGcov = [0,0] # reset

            # dp distrbution of low/high meth CGs
            # TO DO: optimize later
            if not math.isnan(intv.meth[i]):
                if intv.meth[i] <= 0.3:
                    value.nCGLowMeth[dp - 1] += 1
                elif intv.meth[i] >= 0.7:
                    value.nCGHighMeth[dp - 1] += 1
                # nCGMethBin
                kbin = min(19, math.ceil(intv.meth[i] * 20))
                value.nCGMethBin[kbin, min(dp, 30)-1] += 1

        elif intv.cgContext[i] == 'CHG':
            value.nCHG +=1
            value.dpCHG += intv.covMeth[i]
            value.covnCHG[:dp] += 1
            value.meCHG[:dp] += intv.meth[i]

            if intv.base[i] == 'C': 
                value.nCHGW += 1
                value.dpCHGW += intv.covMeth[i]
                value.covnCHGW[:dp] += 1
                value.meCHGW[:dp] += intv.meth[i]
            else:
                value.nCHGC += 1
                value.dpCHGC += intv.covMeth[i]
                value.covnCHGC[:dp] += 1
                value.meCHGC[:dp] += intv.meth[i]
            # nCGMethBin
            if not math.isnan(intv.meth[i]):
                kbin = min(19, math.ceil(intv.meth[i] * 20))
                value.nCHGMethBin[kbin, min(dp, 30)-1] += 1
        elif intv.cgContext[i] == 'CHH':
            value.nCHH +=1
            value.dpCHH += intv.covMeth[i]
            value.covnCHH[:dp] += 1
            value.meCHH[:dp] += intv.meth[i]

            if intv.base[i] == 'C': 
                value.nCHHW += 1
                value.dpCHHW += intv.covMeth[i]
                value.covnCHHW[:dp] += 1
                value.meCHHW[:dp] += intv.meth[i]
            else:
                value.nCHHC += 1
                value.dpCHHC += intv.covMeth[i]
                value.covnCHHC[:dp] += 1
                value.meCHHC[:dp] += intv.meth[i]
            # nCGMethBin
            if not math.isnan(intv.meth[i]):
                kbin = min(19, math.ceil(intv.meth[i] * 20))
                value.nCHHMethBin[kbin, min(dp, 30)-1] += 1

        # error bases for A/T sites (DP>=5)
        elif depth >= 5:
            value.ATdp += depth
            # value.misbase += depth - max(max(intv.covW[:,i]), max(intv.covC[:,i]))
            value.misbase += depth - max(intv.covW[:,i] + intv.covC[:,i])

def update_cgkmer(intv: IntervalCoverage):
    dict_cgkmer = params['dict_cgkmer']
    MAXDEPTH = params['MAXDEPTH']
    
    for i in range(intv.length):
        # if intv.cgContext[i] != 'CG': continue
        if intv.cgkmer[i] == '-': continue

        key = intv.cgkmer[i]
        value = dict_cgkmer.get(key)
        if value is None:
            dict_cgkmer[key] = KmerCov()
            value = dict_cgkmer[key]
            
        # DP = np.arange(start=1, stop=101)
        # MAXDP = 100
        dp = min(MAXDEPTH, intv.covMeth[i])
        value.n += 1
        value.dp += intv.covMeth[i]
        value.cov[:dp] += 1
        value.me[:dp] += intv.meth[i]
        if intv.base[i] == 'C':
            value.nW += 1
            value.dpW += intv.covMeth[i]
            value.covW[:dp] += 1
            value.meW[:dp] += intv.meth[i]
        else:
            value.nC += 1
            value.dpC += intv.covMeth[i]
            value.covC[:dp] += 1
            value.meC[:dp] += intv.meth[i]

def update_lambda(cov: CovLambda, intv: IntervalCoverage):
    for i in range(intv.length):
        if intv.base[i] == 'N': continue

        # bin index
        kbin = intv.bin[i]

        cov.length[kbin] += 1
        dpWC = intv.covWsum[i] + intv.covCsum[i]
        cov.dp[kbin] += dpWC
        cov.dpW[kbin] += intv.covWsum[i]
        cov.dpC[kbin] += intv.covCsum[i]
        if dpWC > 0: cov.cov[kbin] += 1

        # base-wise error base rate
        if dpWC >= 20:
            cov.dp20[kbin] += 1
            cov.misbase[kbin] += (dpWC - max(intv.covW[:,i] + intv.covC[:,i])) / dpWC

        if intv.base[i] == 'C' or intv.base[i] == 'G':
            cov.nC +=1
            if intv.covMeth[i] >= 1:
                cov.covnC += 1
                cov.meC += intv.meth[i]
