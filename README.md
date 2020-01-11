# goVCFtoTable
goVCFtoTable is a tool for parsing a pVCF into a long "melted" tab-delimited format. I'm putting it together
to process pVCF's created by [GLnexus](https://github.com/dnanexus-rnd/GLnexus) for use in a MariaDB database, 
so the tool is a very limited WIP.

Assumes that multiallelic sites are already split.

## Prerequisites
goVCFtoTable uses the bgzf package from biogo for concurrent decompression.
```
$ go get github.com/biogo/biogo/hts/...
```
## Usage
```bash
./goVCFtoTable input.vcf.gz output.tsv numThreads 
```

Comparison with bcftools for chromosome 21 from whole exome sequencing of >83000 people:
```shell script
$ time bcftools query -f '[%CHROM:%POS:%REF:%ALT\t%SAMPLE\t%GT\t%GQ\t%AD\t%FT\n]' -i 'GT="alt"' chr21.vcf.gz > a.txt
real    75m6.131s                                                                                        │
user    74m44.026s                                                                                       │
sys     0m17.645s

$ time ./goVCFtoTable chr21.vcf.gz aaaaaaa.txt 1
real    21m14.092s                                                                                       │
user    27m25.216s                                                                                       │
sys     2m13.799s

$ time ./goVCFtoTable chr21.vcf.gz aaaaaaaaaa.txt 4 
real    14m45.272s                                                                                       │
user    37m10.282s                                                                                       │
sys     2m11.511s
```