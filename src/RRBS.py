

import itertools
from src.utils import *
from src.config import params, data


def compt_RRBS() -> None:
    dict_cgkmer = params['dict_cgkmer']
    dict_rrbs = {
        'ccgg': {'kmer':[], 'n': 0, 'nW': 0, 'nC': 0, 'dp': 0, 'dpW': 0, 'dpC': 0, 'me': 0., 'meW': 0., 'meC': 0., 'cov': 0},
        'non-ccgg': {'kmer':[], 'n': 0, 'nW': 0, 'nC': 0, 'dp': 0, 'dpW': 0, 'dpC': 0, 'me': 0., 'meW': 0., 'meC': 0., 'cov': 0}
    }

    rrbs_kmer = [f'{a}CCGG{b}' for a, b in itertools.product(BASES, BASES)]
    # RRBSKEYS = ['ccgg', 'non-ccgg']

    for key, value in dict_cgkmer.items():
        rrbskey = 'ccgg' if key in rrbs_kmer else 'non-ccgg'
        dict_rrbs[rrbskey]['kmer'].append(key)
        dict_rrbs[rrbskey]['n'] += value.nW + value.nC
        dict_rrbs[rrbskey]['nW'] += value.nW
        dict_rrbs[rrbskey]['nC'] += value.nC
        dict_rrbs[rrbskey]['dp'] += value.dpW + value.dpC
        dict_rrbs[rrbskey]['dpW'] += value.dpW
        dict_rrbs[rrbskey]['dpC'] += value.dpC
        dict_rrbs[rrbskey]['me'] += cumsumrev(value.meW + value.meC)/10000
        dict_rrbs[rrbskey]['meW'] += cumsumrev(value.meW)/10000
        dict_rrbs[rrbskey]['meC'] += cumsumrev(value.meC)/10000
        dict_rrbs[rrbskey]['cov'] += cumsumrev(value.covW + value.covC)

    params['dict_rrbs'] = dict_rrbs
    # params['MAX_DP_RRBS'] = 30
    # maxdp = 30
    return None

def plot_RRBS_CpG_motifs() -> None:
    dict_rrbs = params['dict_rrbs']
    # maxdp = params['MAX_DP_RRBS'] + 1
    maxdp = 30
    img_dir = params['img_dir']

    fig, axs = plt.subplots(2, 2, figsize=(5,4), width_ratios=[1,3],layout='constrained')
    cglabels = ['CCGG', 'non-\nCCGG']
    barcols = [COLS[1], COLS[0]]

    a = dict_rrbs['ccgg']['n']
    b = dict_rrbs['non-ccgg']['n']
    props = [a/(a+b), b/(a+b)]

    data['ccgg_prop'] = fp(props[0])

    bar1 = axs[0,0].bar(cglabels, [a/(a+b), b/(a+b)], color=barcols)
    axs[0,0].set_ylim(top=max(props)+0.1)
    axs[0,0].bar_label(bar1, labels=[f'{e*100:.1f}%' for e in props],
                #  padding=8, color='b', 
                    fontsize=9
                )
    axs[0,0].set_ylabel('proportion\nin genome')

    x = np.arange(1, maxdp)
    axs[0,1].plot(x, dict_rrbs['ccgg']['cov'][1:maxdp]/dict_rrbs['ccgg']['n'], color=COLS[1], label='CCGG')
    axs[0,1].plot(x, dict_rrbs['non-ccgg']['cov'][1:maxdp]/dict_rrbs['non-ccgg']['n'], color=COLS[0], label='non-CCGG')
    axs[0,1].legend()
    axs[0,1].set_ylabel('coverage rate')

    a = dict_rrbs['ccgg']['dp']/dict_rrbs['ccgg']['n']
    b = dict_rrbs['non-ccgg']['dp']/dict_rrbs['non-ccgg']['n']

    data['ccgg_mean_dp'] = ff(a)

    bar2 = axs[1,0].bar(cglabels, [a, b], color=barcols)
    axs[1,0].bar_label(bar2, labels=[f'{e:.1f}' for e in [a, b]],
                #  padding=8, color='b', 
                    fontsize=9
                )
    axs[1,0].set_ylim(top=max([a, b])+1)
    axs[1,0].set_ylabel('mean depth')

    x = np.arange(1, maxdp)
    axs[1,1].plot(x, dict_rrbs['ccgg']['me'][1:maxdp]/dict_rrbs['ccgg']['cov'][1:maxdp], color=COLS[1])
    axs[1,1].plot(x, dict_rrbs['non-ccgg']['me'][1:maxdp]/dict_rrbs['non-ccgg']['cov'][1:maxdp], color=COLS[0])
    axs[1,1].set_xlabel('depth threshold')
    axs[1,1].set_ylabel('mean DNAme level')

    filename = f'{img_dir}/rrbs-kmer-summary'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None
