# ![icon](report/assets/img/favicon-32x32.png) bsDoctor: Quality Diagnosis for Bisulfite-Seq Data

![workflow](report/assets/img/workflow.svg)

## Examples

- Human WGBS sample, https://hippo-yf.github.io/bsDoctor-wgbs-human/
- Human RRBS sample (high double-stranded coverage), https://hippo-yf.github.io/bsDoctor-rrbs-human-high-cov/
- Human RRBS sample (low double-stranded coverage), https://hippo-yf.github.io/bsDoctor-rrbs-human/
- Arabidopsis WGBS sample, https://hippo-yf.github.io/bsDoctor-wgbs-arabidopsis/  

## Features

- define the **diagnostics** of bisulfite-seq data: **bias**, **error**, and **integrity**
- **post-alignment** diagnosis, based-on bam file rather than sequencing fastq file
- **versatile** by a lot of functionalities
- **novel** and significant diagnosis of error, bias, strandness, and CpG motif compared with existing tools
- **flexible** usage to enable/disable each module
- **elaborate** html report for interactive exploration
- pixel (`png`) and **vector (`svg`) figures** for re-edition
- optimized sampling for accuracy-speed **balance**, ~0.5 hours for typical human RRBS sample, ~2 hours for typical human WGBS sample

## Functionalities

- Summary
- Quality
- Bisulfite conversion
- **Base-level error**
- Saturation curve
- Extranuclear DNA: mitochondria/plastid/lambda control DNA
- Genome coverage
- DNAme level
- **DNAme level bias**
- **CpG strandness**
- **CpG motif**
- **RRBS motif**
- Pan-gene DNAme level

## Installation and dependencies

`git clone` or download the repository and run `bsDoctor.py` in any environment with the following dependencies.

- python >= 3.12
- numpy >= 2.1.2 
- pysam >= 0.22.1 
- jinja2 >= 3.1.4
- matplotlib >= 3.9.2
- scipy >= 1.14.1
- tqdm >= 4.66.5

Or use [`uv`](https://docs.astral.sh/uv/) to create virtual environment and install dependencies with the following script. (recommended)  
```shell
cd bsDoctor
uv sync
```

## Usages

- *Reference genome (`-f genome.fa`) must be uncompressed or compressed (`.gz`) with `bgzip`.*

- *Specify the contig names of the mitochondrial (MT), plastid, and (lambda) control DNA if applicable.*

With a downsampled arabidopsis WGBS sample (`examples/example.bam`):

- defaults for all chromosomes and all diagnosis modules:  
`python bsDoctor.py -b examples/example.bam -f examples/genome.fna.gz -g examples/genome.gtf.gz --chr all --mt NC_037304.1 --plastid NC_000932.1 -o output1`

- a subset of chromosomes:  
`python bsDoctor.py -b examples/example.bam -f examples/genome.fna.gz -g examples/genome.gtf.gz --chr NC_003070.9,NC_003071.7  --mt NC_037304.1 --plastid NC_000932.1 -o output2`

- disable pan-gene diagnosis:  
`python bsDoctor.py -b examples/example.bam -f examples/genome.fna.gz --mt NC_037304.1 --plastid NC_000932.1 -o output3 --diag-pangene no`

- disable MT, plastid, and pan-gene diagnosis:  
`python bsDoctor.py -b examples/example.bam -f examples/genome.fna.gz -g examples/genome.gtf.gz --mt NC_037304.1 --plastid NC_000932.1 --diag-mt no --diag-plastid no -o output4`

## Parameters

|**parameter** | **type** | **description**| **defaults** |
|  ----  | ----  | ----  | ----  |
|-b/--bam|str |a .bam file, with **index file** ready|required|
|-f/--fa|str |a .fa[.gz] file of reference genome, uncompressed or compressed with `bgzip`|required|
|-g/--gtf|str |a .gtf[.gz] file|conditional|
|-o/--report-dir|str |report directory|bsDoctor-report|
|--chr|str |nuclear chromosomes to diagnose, separated by comma(,), such as "chr1,chr2,chr3"|all|
|--mt|str |name of mitochondrial DNA|-|
|--plastid|str |name of plastid DNA|-|
|--control-DNA|str |name of spiked-in DNA, usually lambda DNA|-|
|--diag-quality|str |diagnose sequencing and mapping quality or not, true/false, or yes/no|yes|
|--diag-pangene|str |diagnose pangene methylation or not|yes|
|--diag-motif|str |diagnose CpG-motif-related patterns or not|yes|
|--diag-saturation|str |diagnose sequencing saturation or not|yes|
|--diag-mt|str |diagnose reads mapped to mitochondrial DNA or not|yes|
|--diag-plastid|str |diagnose reads mapped to plastid DNA or not|yes|
|--diag-control|str |diagnose reads mapped to spiked-in control DNA (usually lambda DNA) or not|yes|
|--exclude-contigs-start-with|str |exclude contigs, separated by comma(,)|chrUn_|
|--exclude-contigs-end-with|str |exclude contigs, separated by comma(,)|_random,_alt|
|--sampling-prop|float |sampling a proportion of sites of nuclear chromosomes|0.1|
|--bin-size|int |bin size of nuclear chromosomes, overwrite "`--max-chr-nbins`"|0 (use "`--max-chr-nbins`")|
|--max-chr-nbins|int |max number of bins of a nuclear chromosome|1500|
|--bins-control|int |bins of spiked-in control DNA|1000|
|--bins-mt|int |bins of mitochondrial DNA|1000|
|--bins-plastid|int |bins of plastid DNA|1000|
|--base-quality|int |base quality threshold|0|
|--read-quality|int |read mapping quality threshold|0|
|--max-depth|int |depth larger than max depth will be truncated|200|
|--reads-for-quality|int |reads to sample for quality statistics|10000|
|--num-pangene|int |genes to sample for pangene methylation level|1000|
|--max-depth-plot|int |max depth for figures plotted for sites of each specific depth|20|
|--max-depth-motif|int |max depth in CpG-motif diagnosis|50|
|--swap-strand|str |swap read counts on two strands or not|no|
|--ploidy|int |ploidy of the genome|2|
|--save-svg| |save .svg figures or not, save both .png and .svg figures|yes|
|-v/--version| |show program's version number and exit||
|-h/--help|  |show help message and exit||
