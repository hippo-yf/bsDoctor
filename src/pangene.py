import tqdm

from src.utils import *
from src.config import *
from src.coverage import genesGenerator, CovPanGene

def pangene_sampling() -> None:
    gtffile = params['gtffile']
    PANGENE_SAMPLED = params['PANGENE_SAMPLED']

    genes = genesGenerator(gtffile)

    genes_coding = []
    genes_lnc = []
    genes_non = []
    for g in genes:
        if g.end - g.start <= 100: continue # filter short genes
        if g.biotype == 'protein_coding': 
            genes_coding.append(g)
        elif g.biotype == 'lncRNA':
            genes_lnc.append(g)
        else:
            genes_non.append(g)

    genes_coding_sampling = sampleGenes(genes_coding, k=PANGENE_SAMPLED)
    genes_lnc_sampling = sampleGenes(genes_lnc, k=PANGENE_SAMPLED)
    genes_non_sampling = sampleGenes(genes_non, k=PANGENE_SAMPLED)
    
    params['genes_coding_sampling'] = genes_coding_sampling
    params['genes_lnc_sampling'] = genes_lnc_sampling
    params['genes_non_sampling'] = genes_non_sampling

    data['ngene_coding'] = fi(len(genes_coding))
    data['ngene_lnc'] = fi(len(genes_lnc))
    data['ngene_noncoding'] = fi(len(genes_non))
    data['PANGENE_SAMPLED'] = fi(PANGENE_SAMPLED)
    return None

def _plot_pangene_meth(meth: CovPanGene, name: str) -> None:
    gene_breaks = params['gene_breaks']
    img_dir = params['img_dir']

    fig, ax = plt.subplots(figsize=(5,3))
    x = np.arange(gene_breaks*3)
    y = meth.meCG/meth.nCG
    ax.plot(x, y, '-', c=COL_gray)
    ax.plot(x, meth.meCGW/meth.nCGW, '-', c=COLS[1])
    ax.plot(x, meth.meCGC/meth.nCGC, '-', c=COLS[0])
    ax.vlines((49, 99), ax.get_ylim()[0], ax.get_ylim()[1], color= 'gray', linestyles='dashed')
    plt.xticks([0, 49, 74, 99, 149], labels=['-5kbp', 'TSS', 'gene body', 'TTS', '+5kbp'])

    filename = f'{img_dir}/{name}'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

def pangene_compt_plot_meth() -> None:
    genes_coding_sampling = params['genes_coding_sampling']
    genes_lnc_sampling = params['genes_lnc_sampling']
    genes_non_sampling = params['genes_non_sampling']
    bam = params['bam']
    bins = params['gene_breaks']*3

    # coding
    cov_pangene_coding = CovPanGene(bins)
    for i in tqdm.trange(len(genes_coding_sampling)):
        bam.update_pangene(cov_pangene_coding, genes_coding_sampling[i])
    _plot_pangene_meth(cov_pangene_coding, 'pan-gene-meth-coding')

    # lncRNA
    cov_pangene_lnc = CovPanGene(bins)
    for i in tqdm.trange(len(genes_lnc_sampling)):
        bam.update_pangene(cov_pangene_lnc, genes_lnc_sampling[i])
    _plot_pangene_meth(cov_pangene_lnc, 'pan-gene-meth-lncRNA')

    # other noncoding
    cov_pangene_non = CovPanGene(bins)
    for i in tqdm.trange(len(genes_non_sampling)):
        bam.update_pangene(cov_pangene_non, genes_non_sampling[i])
    _plot_pangene_meth(cov_pangene_non, 'pan-gene-meth-other-noncoding')
    return None
