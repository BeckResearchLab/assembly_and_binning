# Metagenomic assembly and binning pipeline for Beck Research Lab

### Generalized protocol
* Trim (see separate protocol repo)
* Initial assembly
  * Use all reads from all libraries
  * PE reads are used as pairs, SE as singletons
* Map raw (untrimmed) reads back to contigs for coverage
* Bin with metabat
  * Successive rounds of binning with parameter searches
  * Choose bins with best precision / recall #s
* Map trimmed reads back to contigs for individual bins
* Extract mapped reads
  * Include pairs, even if unmapped
* Reassemble bins with mapped reads using spades

### Required software
* anaconda3
* bam-readcount
* bedtools
* blast
* bowtie2
* bwa
* checkm
* CONCOCT
* DESMAN
* diamond
* hmmer
* megahit
* metabat
* parallel
* picard
* pplacer
* prodigal
* R
* samtools
* spades
