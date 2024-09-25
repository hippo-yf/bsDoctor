
from numpy import dtype
from src.utils import *
from src.config import params, data

def compt_plot_DNA_composition() -> None:
    binning_covW = params['binning_covW']
    binning_covC = params['binning_covC']
    length = params['length']
    totalBases = params['totalBases']
    img_dir = params['img_dir']

    ## single-stranded coverage rate DP>=1
    single_DP1_coverage = (binning_covW[:,1].sum() + binning_covC[:,1].sum()) / length.sum() / 2

    # unsequenced base prop by Emperical Bayes
    prop_unsequenced =  (binning_covW[:,1].sum() - binning_covW[:,2].sum() + binning_covC[:,1].sum() - binning_covC[:,2].sum()) / totalBases

    # prop of whole genome
    prop_unsequenced = single_DP1_coverage * prop_unsequenced/(1-prop_unsequenced)
    DNA_lost = max(0, 1 - prop_unsequenced - single_DP1_coverage)
    params['DNA_lost'] = DNA_lost
    params['prop_cov_DP1'] = single_DP1_coverage

    ## sequencing saturation 
    saturation_level = 1 - prop_unsequenced

    data['prop_cov_DP1'] = fp(single_DP1_coverage)
    data['prop_unsequenced'] = fp(prop_unsequenced)
    data['prop_asymptotic'] = fp(single_DP1_coverage + prop_unsequenced)
    data['prop_lost'] = fp(DNA_lost)
    data['saturation'] = fp(saturation_level)

    # pie plot of DNA
    fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))

    DNA_composition = [
        "sequenced (depth ≥ 1)",
        "undegraded but not sequenced",
        "degraded"
        ]

    props = [single_DP1_coverage, prop_unsequenced, DNA_lost]
    wedges, texts, autotexts = ax.pie(
        props,
        # labels=DNA_composition, 
        autopct=lambda x:fp(x/100)+'%', 
        textprops=dict(color="w"),
        colors=COLS[:3]
        )

    ax.legend(wedges, DNA_composition,
            title="",
            loc="center left",
            bbox_to_anchor=(1, 0, 0.5, 1))
    plt.setp(autotexts, size=8, weight="bold")

    filename = f'{img_dir}/pie-DNA-lost'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

def dilute_bin(cov: NDArray) -> Tuple[NDArray, NDArray]:
    """
    dilution of genomic bin
    cov: cummulative counts 
    """
    L = len(cov)
    DPS = [0,1,3,5,10]
    maxdp = max(max(DPS)+1, L - np.argmax(cov[::-1] > 0))
    cov = depthDiff(cov[:maxdp])

    dps = np.arange(maxdp)
    tbases = dps @ cov
    breaks = 50
    dilution = np.zeros((breaks+1, len(DPS)), dtype=int32) # 30bins X dp=[0,1,3,5,10]
    dilution[0,0] = cov.sum()
    residue = tbases
    bases_to_sample = np.array(tbases*(np.arange(0, breaks+1)/breaks), dtype=int64)
    k = breaks

    for _ in range(tbases):
        if residue <= bases_to_sample[k]:
            dilution[k,:] = cumsumrev(cov)[DPS]
            k -= 1
        if k < 1: break
        weights = cov * dps
        weights = weights/weights.sum()
        i = np.where(GEN.multinomial(1, weights))[0]
        cov[i] -= 1
        cov[i-1] += 1
        residue -=1 
    if k > 0:
        dilution[1:(k+1),0] = dilution[0,0]
    return (bases_to_sample, dilution)

def dilute_strand(cov: NDArray, nbins: int = 100) -> Tuple[NDArray, NDArray]:
    """
    dilution of bins of singe strand
    """
    nbases = np.zeros((51,), dtype=float64)
    dilu = np.zeros((51, 5), dtype=int64)

    shape = cov.shape
    if shape[0] > nbins *1.1:
        kbins = GEN.choice(shape[0], size = nbins, replace=False)
    else:
        kbins = range(shape[0])
    for i in kbins:
        bases_to_sample, dilution = dilute_bin(cov[i,:])
        nbases += bases_to_sample
        dilu += dilution
    return (nbases, dilu)

def dilute_genome(covW: NDArray, covC: NDArray, nbins: int = 100) -> Tuple[NDArray, NDArray]:
    """
    dilution of two strands
    """
    nbasesW, diluW = dilute_strand(covW, nbins=nbins)
    nbasesC, diluC = dilute_strand(covC, nbins=nbins)
    nbases = nbasesW + nbasesC
    dilu = diluW + diluC
    return (nbases, dilu)

def plot_saturation_curve() -> None:
    bam = params['bam']
    DNA_lost = params['DNA_lost']
    binning_covW = params['binning_covW']
    binning_covC = params['binning_covC']
    prop_cov_DP1 = params['prop_cov_DP1']
    img_dir = params['img_dir']
    mean_read_length = params['mean_read_Mlen']

    library_bases = mean_read_length * bam.mapped

    dps = [1,3,5,10]
    cols = ['#a0d0f8', '#40a0f2','#0d6dbf', '#07375f']

    (nbases, dilu) = dilute_genome(binning_covW, binning_covC, nbins=100)

    fig, ax = plt.subplots(figsize=(5, 3))
    prop_dp1_subset = dilu[50,1]/dilu[50,0]

    x = np.array(nbases/nbases[-1]*library_bases, dtype=int64)
    x, suffix = getSuffix(x)
    for i in range(len(dps)):
        y = dilu[:,i+1]/dilu[:,0] * 100
        # rescale to whole-genole covrate
        y = y*prop_cov_DP1/prop_dp1_subset 
        ax.plot(x, y, '-', c=cols[i], alpha=1, label='DP>='+str(dps[i]))
    abline((1-DNA_lost)*100, 0, color=COL_gray)
    ax.axis()
    linewidth = 0.001
    ax.text(0, (1-DNA_lost-linewidth)*100, f'{100-DNA_lost*100:.2f}%', verticalalignment='top')
    ax.legend(title='single-stranded')
    plt.xlabel(f'sequenced bases ({suffix} bases)')
    # plt.xlabel('single-stranded sequencing depth')
    plt.ylabel('single-stranded coverage (%)')

    filename = f'{img_dir}/saturation-curve'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    