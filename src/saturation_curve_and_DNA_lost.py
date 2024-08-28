
from src.utils import *
from src.config import params, data

def compt_plot_DNA_composition() -> None:
    binning_covW = params['binning_covW']
    binning_covC = params['binning_covC']
    length = params['length']
    totalBases = params['totalBases']
    img_dir = params['img_dir']

    ## coverage rate single-strand DP>=1
    single_DP1_coverage = (binning_covW[:,0].sum() + binning_covC[:,0].sum()) / length.sum() / 2

    # unsequenced base prop by Emperical Bayes
    prop_unsequenced =  (binning_covW[:,0].sum() - binning_covW[:,1].sum() + binning_covC[:,0].sum() - binning_covC[:,1].sum()) / totalBases

    # prop of whole genome
    prop_unsequenced = single_DP1_coverage * prop_unsequenced/(1-prop_unsequenced)
    DNA_lost = max(0, 1 - prop_unsequenced - single_DP1_coverage)
    params['DNA_lost'] = DNA_lost

    ## sequencing saturation 
    saturaion_level = 1 - prop_unsequenced

    data['prop_cov_DP1'] = fp(single_DP1_coverage)
    data['prop_unsequenced'] = fp(prop_unsequenced)
    data['prop_asymptotic'] = fp(single_DP1_coverage + prop_unsequenced)
    data['prop_lost'] = fp(DNA_lost)
    data['saturaion'] = fp(saturaion_level)

    # pie plot of DNA
    fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))

    DNA_composition = [
        "sequenced (depth ≥ 1)",
        "undegraded but not sequenced",
        "degenerated"
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
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

def _saturation_curve_binning(sampled_bases: int, cov: NDArray) -> NDArray:

    maxdp = len(cov) - np.argmax(cov[::-1] > 0)
    depths = np.repeat(np.arange(1, maxdp+1), cov[:maxdp])

    if sampled_bases >= np.sum(depths):
        return np.array(cov[ [0,2,4,9] ])
    sampling = GEN.multivariate_hypergeometric(depths, sampled_bases)

    return np.array([np.sum(sampling>=1), np.sum(sampling>=3), np.sum(sampling>=5), np.sum(sampling>=10)])

def _saturation_curve(breaks=30) -> NDArray:
    """
    single-stranded saturation curve
    """
    binning_depthW = params['binning_depthW']
    binning_depthC = params['binning_depthC']
    binning_covW = params['binning_covW']
    binning_covC = params['binning_covC']
    length = params['length']
    totalBasesW = params['totalBasesW']
    totalBasesC = params['totalBasesC']
    totalDepthW = params['totalDepthW']
    totalDepthC = params['totalDepthC']
    
    genome_coverages = np.zeros(shape=(breaks+1, 4))
    genome_coverages[breaks,:] = np.array([
        binning_covW[:,0].sum() + binning_covC[:,0].sum(),
        binning_covW[:,2].sum() + binning_covC[:,2].sum(),
        binning_covW[:,4].sum() + binning_covC[:,4].sum(),
        binning_covW[:,9].sum() + binning_covC[:,9].sum()]
        )/length.sum()/2
    
    # sampled_bases = np.logspace(0, np.log10(totalDepth.sum()), breaks+1, dtype=np.int64) -1
    sampled_basesW = np.linspace(0, totalBasesW, breaks+1, dtype=np.int64)
    sampled_basesC = np.linspace(0, totalBasesC, breaks+1, dtype=np.int64)
    for i in np.arange(1, breaks):

        # Watson strand
        bin_basesW = GEN.multinomial(sampled_basesW[i], totalDepthW/totalBasesW, size=1)[0,:]
        # DP >= 1,3,5,10
        coveraged_basesW = [
            _saturation_curve_binning(bin_basesW[j], binning_depthW[j,:])
            for j in range(len(length))
            ]

        # Crick strand
        bin_basesC = GEN.multinomial(sampled_basesC[i], totalDepthC/totalBasesC, size=1)[0,:]
        # bin_bases = np.fmin(bin_bases, totalDepth)

        # DP >= 1,3,5,10
        coveraged_basesC = [
            _saturation_curve_binning(bin_basesC[j], binning_depthC[j,:])
            for j in range(len(length))
            ]
        
        genome_coverages[i,:] = (np.asarray(coveraged_basesW).sum(axis=0) + np.asarray(coveraged_basesC).sum(axis=0))/length.sum()/2

    # add mean sequencing depth
    return np.hstack([
        genome_coverages, 
        (sampled_basesW+sampled_basesC).reshape((-1,1))/length.sum()/2
        ])

def plot_saturation_curve() -> None:
    DNA_lost = params['DNA_lost']
    img_dir = params['img_dir']

    sc = _saturation_curve(breaks=30)

    dp = [1,3,5,10]
    cols = ['#a0d0f8', '#40a0f2','#0d6dbf', '#07375f']

    fig, ax = plt.subplots(figsize=(5, 3))
    size = np.shape(sc)
    x = sc[:,-1]
    for i in range(size[1]-1):
        ax.plot(x, sc[:,i], '-', c=cols[i], alpha=1, label='DP>='+str(dp[i]))
    # ax.plot(np.linspace(0, sc[-1, -1], 10), np.repeat(1-DNA_lost, 10), '--')
    abline(1-DNA_lost, 0, color=COL_gray)
    linewidth = 0.001
    ax.text(0, 1-DNA_lost-linewidth, f'{100-DNA_lost*100:.2f}%', verticalalignment='top')
    ax.legend(title='single-stranded')
    plt.xlabel('single-stranded sequencing depth')
    plt.ylabel('single-stranded genomic coverage')

    filename = f'{img_dir}/saturation-curve'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None
    