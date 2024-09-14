

from src.config import params, data
from src.utils import *
from src.coverage import randomGenomicRegion, Quality

#### without sampling

def compt_plot_DNA_content() -> None:
    bam = params['bam']
    fa = params['fa']
    # reference_length = params['reference_length']
    reference_length = bam.get_reference_length
    chr_MT = params['chr_MT']
    chr_lambda = params['chr_lambda']
    chr_plastid = params['chr_plastid']
    img_dir = params['img_dir']
    ploidy = params['ploidy']

    # reads of nuclear, mt, plastid, and lambda
    nmapped = {'nuclear': 0, 'mt': 0, 'lambda': 0, 'plastid': 0}
    length = {'nuclear': 0, 'mt': 0, 'lambda': 0, 'plastid': 0} # chr length
    # nmapped = {'nuclear': 0}
    # if chr_MT != '-': nmapped['mt'] = 0
    # if chr_lambda != '-': nmapped['lambda'] = 0
    # if chr_plastid != '-': nmapped['plastid'] = 0

    # reads of each chr
    nmapped_chrs = dict()

    for c in bam.get_index_statistics():
        nmapped_chrs[c.contig] = c.mapped
        if chr_MT != '-' and c.contig == chr_MT: 
            nmapped['mt'] += c.mapped
            length['mt'] += reference_length(c.contig)
        elif chr_lambda != '-' and c.contig == chr_lambda: 
            nmapped['lambda'] += c.mapped
            length['lambda'] += reference_length(c.contig)
        elif chr_plastid != '-' and c.contig == chr_plastid: 
            nmapped['plastid'] += c.mapped
            length['plastid'] += reference_length(c.contig)
        else: 
            nmapped['nuclear'] += c.mapped
            length['nuclear'] += reference_length(c.contig)
        
    params['nmapped_chrs'] = nmapped_chrs

    # read length from a samll batch of reads
    genomeGen = randomGenomicRegion(fa, bam, size=10)
    quality = Quality()
    isampled = 0
    while quality.nreads < 10**3 and isampled < 10_000:
        region = next(genomeGen)
        bam.sample_reads(quality=quality, gi=region)
        isampled += 1
    quality.regularization()

    mean_read_len = np.mean(quality.read_length)

    ## nreads and prop in report table
    _nreads = dict()
    _nreads_prop = dict()
    _mean_dp = dict()
    _copy_num = dict()
    for key, value in nmapped.items():
        _nreads[key] = fi(value)
        _nreads_prop[key] = fp(value/bam.mapped)
        _mean_dp[key] = fi(int(mean_read_len*value/length[key])) if length[key]>0 else "0"
        _copy_num[key] = fi(int((nmapped[key]/length[key]) / (nmapped['nuclear']/length['nuclear']) * ploidy)) if length[key]>0 else "0"
    data['nreads'] = _nreads
    data['nreads_prop'] = _nreads_prop
    data['mean_dp_dict'] = _mean_dp
    data['copy_num'] = _copy_num


    lambda_DNA_content_vs_nuclear = nandivide(nmapped['lambda'], nmapped['nuclear'])
    lambda_DNA_content_vs_MT = nandivide(nmapped['lambda'], nmapped['mt'])
    lambda_DNA_content_vs_plastid = nandivide(nmapped['lambda'], nmapped['plastid'])
    
    data['lambda_DNA_content_vs_nuclear'] = fp(lambda_DNA_content_vs_nuclear)
    data['lambda_DNA_content_vs_MT'] = fp(lambda_DNA_content_vs_MT)
    data['lambda_DNA_content_vs_plastid'] = fp(lambda_DNA_content_vs_plastid)

    #### reads mapped to each chr
    # for interactive scatter plot
    chrs_fa = fa.references
    nreads = []
    lens = []
    chrs = []
    for c in bam.get_index_statistics():
        cont = c.contig
        if contig_should_be_included(cont) and cont in chrs_fa:
        # if cont in chrs_fa:
            chrs.append(cont)
            nreads.append(c.mapped)
            lens.append(reference_length(cont))
    dict_reads = {'chrs': chrs, 'lens': lens, 'nreads': nreads}
    data['dict_reads'] = dict_reads

    #### pie plot
    
    labels = {'nuclear': 'nuclear', 'mt': 'MT', 'plastid': 'plastid', 'lambda': 'lambda/spiked-in'}

    nreads = []
    labs = []
    for key, value in nmapped.items():
        if value > 0:
            labs.append(labels[key])
            nreads.append(value)

    # props = np.array([nmapped['nuclear'], nmapped['mt'], nmapped['lambda']]) 
    props = np.array(nreads)
    props = props / np.sum(props)

    fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
    wedges, texts, autotexts = ax.pie(
        props,
        autopct=lambda x:fp(x/100)+'%', 
        textprops=dict(color="w"),
        colors=COLS[:3]
        )
    ax.legend(wedges, labs,
            title="",
            loc="center left",
            bbox_to_anchor=(1, 0, 0.5, 1))
    plt.setp(autotexts, size=8, weight="bold")

    filename = f'{img_dir}/pie-DNA-content'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

