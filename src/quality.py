
from src.config import params, data
from src.utils import *
from src.coverage import randomGenomicRegion
# from src.sampling import *

## summary
def compt_quality() -> None:
    bam = params['bam']
    fa = params['fa']
    reads_sample_window = params['reads_sample_window']
    reads_to_sample = params['reads_to_sample']
    quality = params['quality']
    
    # gi = GenomicInterval('chr1', 2000000, 1000_000, 1001_000)
    # bam.sample_reads(quality, gi, params)
    # quality.regularization()

    genomeGen = randomGenomicRegion(fa, bam, size=reads_sample_window)
    isampled = 0 # sampled reads
    while quality.nreads < reads_to_sample and isampled < 1_000_000:
        region = next(genomeGen)
        bam.sample_reads(quality=quality, gi=region)
        isampled += 1
    quality.regularization()

    params['mean_read_length'] = np.mean(quality.read_length)
    data['quality_nreads'] = fi(quality.nreads)
    data['quality_nbases'] = fi(len(quality.base_quality))
    data['quality_nloci'] = fi(isampled)
    return None


########################
#### figures
########################

def plot_read_length() -> None:
    quality = params['quality']

    try:
        fig, ax = plt.subplots(figsize=(4,2))
        plt.hist(quality.read_length, bins=21, density=True, color=COLS[0])
        plt.xlabel('read length')
        plt.ylabel('density')

        # filename = r'img/dist-read-length'
        filename = params['img_dir'] + r'/hist-read-length'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    except:
        pass
    return None

def plot_base_quality() -> None:
    quality = params['quality']

    try:
        fig, ax = plt.subplots(figsize=(4,2))
        plt.hist(quality.base_quality, bins=21, density=True, color=COLS[0])
        plt.xlabel('base quality')
        plt.ylabel('density')

        filename = params['img_dir'] + r'/hist-base-quality'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    except:
        pass
    return None

def plot_read_map_quality() -> None:
    quality = params['quality']

    try:
        fig, ax = plt.subplots(figsize=(4,2))
        plt.hist(quality.map_quality, bins=21, density=True, color=COLS[0])
        plt.xlabel('read mapping quality')
        plt.ylabel('density')

        filename = params['img_dir'] + r'/hist-read-map-quality'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    except:
        pass
    return None

def plot_bar_base_cigar() -> None:
    quality = params['quality']
    prop = quality.ncigars/quality.ncigars.sum()
    x = range(len(prop))
    CIGARS = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X', 'B']

    try:
        fig, ax = plt.subplots(figsize=(4,2))
        plt.bar(x, prop, color=COLS[0])
        plt.xticks(x, labels=CIGARS)
        plt.xlabel('base mapping CIGAR')
        plt.ylabel('proportion')

        filename = params['img_dir'] + r'/bar-base-CIGAR'
        plt.savefig(filename+'.png', transparent=True, dpi=300, bbox_inches='tight')
        if params['save_svg']: plt.savefig(filename+'.svg', transparent=True, bbox_inches='tight')
        plt.close()
    except:
        pass
    return None
    