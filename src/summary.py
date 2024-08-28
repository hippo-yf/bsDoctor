

from src.config import params, data
from src.utils import *

#### without sampling

def compt_plot_DNA_content() -> None:
    bam = params['bam']
    chr_MT = params['chr_MT']
    chr_lambda = params['chr_lambda']
    img_dir = params['img_dir']

    nmapped = {'nuclear': 0, 'mt': 0, 'lambda': 0}
    for c in bam.get_index_statistics():
        if c.contig == chr_MT: nmapped['mt'] += c.mapped
        elif c.contig == chr_lambda: nmapped['lambda'] += c.mapped
        else: nmapped['nuclear'] += c.mapped

    lambda_DNA_content_vs_nuclear = nmapped['lambda']/nmapped['nuclear']
    lambda_DNA_content_vs_MT = nmapped['lambda']/nmapped['mt']

    data['lambda_DNA_content_vs_nuclear'] = fp(lambda_DNA_content_vs_nuclear)
    data['lambda_DNA_content_vs_MT'] = fp(lambda_DNA_content_vs_MT)

    #### pie plot
    labels = [
        "nucler",
        "MT",
        "lambda/spiked-in"
        ]
    props = np.array([nmapped['nuclear'], nmapped['mt'], nmapped['lambda']]) 
    props = props / np.sum(props)

    fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
    wedges, texts, autotexts = ax.pie(
        props,
        autopct=lambda x:fp(x/100)+'%', 
        textprops=dict(color="w"),
        colors=COLS[:3]
        )
    ax.legend(wedges, labels,
            title="",
            loc="center left",
            bbox_to_anchor=(1, 0, 0.5, 1))
    plt.setp(autotexts, size=8, weight="bold")

    filename = f'{img_dir}/pie-DNA-content'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None
