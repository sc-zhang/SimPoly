## Introduction

SimPoly is a simple tool to simulate polyploid geome from mono genome.

## Dependencies

### Software

- python >= 3.7

### Python packages

- numpy

## Installation

```bash
cd /path/to/install
git clone https://github.com/sc-zhang/SimPoly.git
pip install -r requirements.txt
chmod +x SimPoly/sim_poly.py
# Optional
echo 'export PATH=/path/to/install/SimPoly:$PATH' >> ~/.bash_profile
source ~/.bash_profile
```

## Usage

- Example

```bash
sim_poly.py -g genome.fa -f genome.gff3 -c config.conf -o sim
```

- Details

```bash
usage: sim_poly.py [-h] -g GENOME -f GFF3 -c CONFIG -o OUT

options:
  -h, --help            show this help message and exit
  -g GENOME, --genome GENOME
                        input genome file, gz supported
  -f GFF3, --gff3 GFF3  input gff3 file, gz supported
  -c CONFIG, --config CONFIG
                        input config file
  -o OUT, --out OUT     output directory
```

- Details of config file (example config file could be found at root directory of this tool)

> 1. ploidy is the target ploidy need simulate, for each haplotype, the ratio of SNP,insertion,deletion shoud be set
     > seperated by comma.
> 2. the max length of insertion/deletion fragments are use to generate insert fragments which length
     > are in range of [1, max_length_of_insertion], and delete regions in range of [1, max_length_of_deletion]
> 3. the structure is used to define how to simulate haplotypes
> ```text
> [parameters]
> ploidy=4
> snp=1e-4,1e-4,1e-4,1e-4
> insertion=5e-5,5e-5,5e-5,5e-5
> deletion=5e-5,5e-5,5e-5,5e-5
> # max length of insertion/deletion fragments
> insertion_length=1e3,1e3,1e3,1e3
> deletion_length=1e3,1e3,1e3,1e3
> # structure string is constructed by numbers
> # 0 means reference
> # 1 means hap1
> # 2 means hap2
> # ...
> # 0,0,0,0 means all 4 haploids are simulated from reference sequence
> # 0,1,2,3 means hap1 is simulated from reference, hap2 is simulated from hap2, ...
> # 0,1,0,2 means hap1 and hap3 are simulated from reference, hap2 is simulated from hap1, hap4 is simulated from hap3
> # ...
> structure=0,1,2,3
> 
> ```

## Results

1. snp*.txt, ins*.txt, del*.txt record the detail information of all variants.
2. sim.fasta is the final genome file.
3. sim.gff3 is the final annotation file.