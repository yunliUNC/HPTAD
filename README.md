# HPTAD user manual
## What is HPTAD?
HPTAD is a method for calling topologically associating domains (TADs) from HiChIP or PLAC-seq datasets. 
 
Two of the required inputs for HPTAD can be generated using [<em>MAPS</em>](https://github.com/HuMingLab/MAPS) (detailed in [this paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006982)). This preprocessing pipeline, called <em>feather</em>, converts aligned, sorted, and merged paired-end reads to long-range and short-range .bed/.bedpe files. For questions regarding HPTAD please email Ming Hu (hum@ccf.org) or Yun Li (yunli@med.unc.edu).

The requirements and details for running this software are provided below:

The HPTAD R program was tested and run on Linux and requires several readily available packages. We list the versions used in testing and development with the caveat that (slightly) older and newer versions should work, but caveat utilitor.

**R 3.6.0**
* argparse 2.0.3
* data.table 1.14.2
* dplyr 1.0.7
* MASS 7.3.51.4
* VGAM 1.1.3
* zoo 1.8.5

## Inputs
1. .long.intra.bedpe file: obtained from [<em>MAPS</em>](https://github.com/HuMingLab/MAPS) as referenced above
2. .shrt.vip.bed file: obtained from [<em>MAPS</em>](https://github.com/HuMingLab/MAPS) as referenced above
3. Genomic features file: can be downloaded from [Genomic features](http://enhancer.sdsc.edu/yunjiang/resources/genomic_features/)
4. 1-D ChIP peaks: can be obtained by running <em>MACS2</em> on corresponding ChIP-seq data, or on short-range reads
5. Filter file (optional):
* This is used if you want to exclude genomic regions from MAPS analysis. Reads mapping to those regions will be ignored. Set to “None” if not filtering any regions. A filter file is a tab delimited table containing three columns (no header) representing chromosome, the start position, and end position of any loci you wish to exclude.
```
chr10	22142530 22142880
chr10	22142830 22143070
chr10	35110060 35110270
chr10	58223870 58224100
chr11	39148660 39148860
```

## Usage
Running the script HPTAD.R involves supplying arguments via the command line, most indicating the location of the inputs listed above.
* -i, --indir: directory containing both long-range and short-range .bed/.bedpe files for desired chromosomes (required)
* -o, --outdir: directory where output is to be written (required)
* -p, --prefix: name of dataset: should be identical to prefix of .bed/.bedpe files (required)
* -C, --chip: full path to ChIP peak file (required)
* -c, --chromosomes: comma separated list of chromosomes to be analyzed (1,2,3 e.g.) (required)
* -f, --features: full path to features file (required)
* -x, --filter: full path to filter file, or 'None' (optional: default = 'None')
* -b, --binsize: bin size for analysis in base pair length (required)
* -u, --upperlimit: upper limit of distance between bins to be analyzed (optional: default = 2Mb)

## Example
The files required to run HPTAD on an example dataset are provided in the example folder. Sample command line submission:
```
Rscript HPTAD.R \
  -i example \
  -o test \
  -p foo \
  -C example/ChIP/foo.narrowPeak \
  -c 19 \
  -f example/features/foo_F_GC_M_MboI_40Kb_el.txt \
  -x example/filter/foo_filterlist.bed \
  -b 40000
```
Output is provided as a tab separated value files containing chromosome, start, and end position of TAD regions for each chromosome analyzed.
