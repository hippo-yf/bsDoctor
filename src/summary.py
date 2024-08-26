

from src.config import params, data
from src.utils import *
# from src.coverage import *
# from src.sampling import *

def compt_DNA_content():
    bam = params['bam']
    chr_MT = params['chr_MT']
    chr_lambda = params['chr_lambda']

    nmapped = {'nuclear': 0, 'mt': 0, 'lambda': 0}
    for c in bam.get_index_statistics():
        if c.contig == chr_MT: nmapped['mt'] += c.mapped
        elif c.contig == chr_lambda: nmapped['lambda'] += c.mapped
        else: nmapped['nuclear'] += c.mapped

    lambda_DNA_content_vs_nuclear = nmapped['lambda']/nmapped['nuclear']
    lambda_DNA_content_vs_MT = nmapped['lambda']/nmapped['mt']

    data['lambda_DNA_content_vs_nuclear'] = fp(lambda_DNA_content_vs_nuclear)
    data['lambda_DNA_content_vs_MT'] = fp(lambda_DNA_content_vs_MT)

