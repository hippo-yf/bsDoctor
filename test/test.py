
from argparse import ArgumentParser

class MyArgumentParser(ArgumentParser):
    def __init__(self, **args) -> None:
        super().__init__(**args)
        self.add_argument('-i', '--atcg-file', dest='infile', help='an input .ATCGmap[.gz] file, default: read from stdin', type=str, required=False, default="-")


desc = 'genotype/SNP calling from bisulite-sequencing data'

parser = MyArgumentParser(description=desc)
options = parser.parse_args()
print(options.__dict__)

params = {}

for (key, value) in options.__dict__.items():
    params[key] = value

print(params)
