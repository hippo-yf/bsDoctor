
import numpy as np
from src.config import params, data
from src.utils import *

from importlib import reload
import src.utils
reload(src.utils)
relative_freq_CpG_motif = src.utils.relative_freq_CpG_motif


def compt_CpG_motif() -> None:
    dict_cgkmer =  params['dict_cgkmer']
    maxdp = params['MAX_DP_CG_MOTIF']

    kmer = list()
    kmer_n = list()
    kmer_nW = list()
    kmer_nC = list()
    kmer_dp = list()
    kmer_dpW = list()
    kmer_dpC = list()
    kmer_meth = list()
    kmer_methW = list()
    kmer_methC = list()
    kmer_cov = list()

    for key, value in dict_cgkmer.items():
        kmer.append(key)
        kmer_n.append(value.n)
        kmer_nW.append(value.nW)
        kmer_nC.append(value.nC)
        kmer_dp.append(value.dp / value.n)
        kmer_dpW.append(value.dpW / value.nW)
        kmer_dpC.append(value.dpC / value.nC)
        kmer_meth.append(nandivide(value.me, value.cov))
        kmer_methW.append(nandivide(value.meW, value.covW))
        kmer_methC.append(nandivide(value.meC, value.covC))
        kmer_cov.append(value.cov / value.n)

    # df = pd.DataFrame({'kmer':kmer, 'n':n, 'nW':nW, 'nC':nC, 'dp':dp, 'dpW':dpW, 'dpC':dpC})
    dict_kmer = {
        'kmer':kmer, 
        'n': np.asarray(kmer_n), 
        'nW': np.asarray(kmer_nW), 
        'nC': np.asarray(kmer_nC), 
        'dp': np.asarray(kmer_dp), 
        'dpW': np.asarray(kmer_dpW), 
        'dpC': np.asarray(kmer_dpC),
        'me': np.asarray(kmer_meth),
        'meW': np.asarray(kmer_methW),
        'meC': np.asarray(kmer_methC),
        'cov': np.asarray(kmer_cov)
        }
    params['dict_kmer'] = dict_kmer # a summary

    #### data for interactive plots
    
    # cgkmer_2strands = {
    #     'nW': list(signif(relative_freq_CpG_motif((dict_kmer['nW'])))),
    #     'nC': list(signif(relative_freq_CpG_motif((dict_kmer['nC'])))),
    #     'n': list(signif(relative_freq_CpG_motif((dict_kmer['n'])))),
    #     'cgkmer': dict_kmer['kmer'],
    #     'mean_dp': list(signif((dict_kmer['dp']), 2)),
    #     'meWdp1': list(signif((dict_kmer['meW'][:,0]))),
    #     'meCdp1': list(signif((dict_kmer['meC'][:,0]))),
    #     'medp1': list(signif((dict_kmer['me'][:,0]))),
    #     'meW': list(signif((dict_kmer['meW'][:,0]))),
    #     'meC': list(signif((dict_kmer['meC'][:,0]))),
    #     'dp': list(range(1, maxdp + 1)),
    #     'cov': [list(signif((x))) for x in dict_kmer['cov'][:,:maxdp]],
    #     'me': [list(signif((x))) for x in dict_kmer['me'][:,:maxdp]]
    # }
    cgkmer_2strands = {
        'nW': [float(signif(x)) for x in relative_freq_CpG_motif(dict_kmer['nW'])],
        'nC': [float(signif(x))for x in relative_freq_CpG_motif(dict_kmer['nC'])],
        'n': [float(signif(x)) for x in relative_freq_CpG_motif(dict_kmer['n'])],
        'cgkmer': dict_kmer['kmer'],
        'mean_dp': [float(signif(x)) for x in dict_kmer['dp']],
        'meWdp1': [float(signif(x)) for x in dict_kmer['meW'][:,0]],
        'meCdp1': [float(signif(x)) for x in dict_kmer['meC'][:,0]],
        'medp1': [float(signif(x)) for x in dict_kmer['me'][:,0]],
        # 'meW': list(signif((dict_kmer['meW'][:,0]))),
        # 'meC': list(signif((dict_kmer['meC'][:,0]))),
        # 'dp': list(range(1, maxdp + 1)),
        # 'cov': [list(signif((x))) for x in dict_kmer['cov'][:,:maxdp]],
        'cov': [[float(signif(y)) for y in x if not np.isnan(y)] for x in dict_kmer['cov'][:,:maxdp]],
        'me': [[float(signif(y)) for y in x if not np.isnan(y)] for x in dict_kmer['me'][:,:maxdp]],
    }
    data['cgkmer_2strands'] = cgkmer_2strands
    return None

def plot_hist_CpG_motif_freq() -> None:
    dict_kmer = params['dict_kmer']
    img_dir = params['img_dir']

    fig, ax = plt.subplots(figsize=(5,3))
    ax.hist(relative_freq_CpG_motif(dict_kmer['n']), bins=21, color=COLS[0])
    ax.vlines(1, *ax.get_ylim(), linestyles='dashed', color=COL_gray)
    plt.xlabel('relative freq (fold over 1/256)')
    plt.ylabel('#CpG-motif')

    filename = f'{img_dir}/kmer-frequency-dist'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

def plot_CpG_motif_freq_watson_vs_crick() -> None:
    dict_kmer = params['dict_kmer']
    img_dir = params['img_dir']

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot(relative_freq_CpG_motif(dict_kmer['nW']), relative_freq_CpG_motif(dict_kmer['nC']), 'o', c=COLS[0], markersize=2)
    abline(0, 1, color=COL_gray)
    plt.xlabel('relative freq on Watson strand')
    plt.ylabel('relative freq on Crick strand')

    filename = f'img/{img_dir}-stranded-frequency'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

def plot_CpG_motif_me_watson_vs_crick() -> None:
    dict_kmer = params['dict_kmer']
    img_dir = params['img_dir']

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot(dict_kmer['meW'][:,0], dict_kmer['meC'][:,0], 'o', c=COLS[0], markersize=2)
    abline(0, 1, color=COL_gray)
    plt.xlabel('mean DNAme level on Watson strand')
    plt.ylabel('mean DNAme level on Crick strand')

    filename = f'{img_dir}/kmer-stranded-me-dp1'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

def plot_CpG_motif_depth_vs_freq() -> None:
    dict_kmer = params['dict_kmer']
    img_dir = params['img_dir']

    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot(relative_freq_CpG_motif(dict_kmer['n']), dict_kmer['dp'], 'o', c=COLS[0], markersize=2)
    plt.xlabel("relative freq on two strands")
    plt.ylabel("mean read depth")
    plt.vlines(1, ax.get_ylim()[0], ax.get_ylim()[1], colors=COL_gray, linestyles='dashed')

    filename = f'{img_dir}/kmer-mean-depth-vs-freq'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

def plot_CpG_motif_me_vs_freq() -> None:
    dict_kmer = params['dict_kmer']
    img_dir = params['img_dir']
        
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.plot(relative_freq_CpG_motif(dict_kmer['n']), dict_kmer['me'][:,0], 'o', c=COLS[0], markersize=2)

    plt.xlabel("relative freq on two strands")
    plt.ylabel("mean DNAme level")
    plt.vlines(1, ax.get_ylim()[0], ax.get_ylim()[1], colors=COL_gray, linestyles='dashed')

    filename = f'{img_dir}/kmer-meth-vs-frequency'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

def plot_CpG_motif_covrate_vs_depth() -> None:
    dict_kmer = params['dict_kmer']
    img_dir = params['img_dir']
    maxdp = params['MAX_DP_CG_MOTIF']

    fig, ax = plt.subplots(figsize=(5, 3))
    dp = np.arange(start=1, stop=maxdp+1)
    for i in range(np.shape(dict_kmer["cov"])[0]):
        ax.plot(dp, dict_kmer['cov'][i,:maxdp], '--', c=COLS[0], alpha=0.15)
    plt.xlabel("depth threshold")
    plt.ylabel("coverage rate")

    filename = f'{img_dir}/kmer-coverage-vs-depth'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

def plot_CpG_motif_me_vs_depth() -> None:
    dict_kmer = params['dict_kmer']
    img_dir = params['img_dir']
    maxdp = params['MAX_DP_CG_MOTIF']
        
    fig, ax = plt.subplots(figsize=(5, 3))
    dp = np.arange(start=1, stop=maxdp+1)
    for i in range(np.shape(dict_kmer['cov'])[0]):
        ax.plot(dp, dict_kmer['me'][i,:maxdp], '--', c=COLS[0], alpha=0.15)
    plt.xlabel("depth threshold")
    plt.ylabel("mean DNAme level")

    filename = f'{img_dir}/kmer-meth-vs-depth'
    plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
    plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
    plt.close()
    return None

