echo "mapping reads with $@"

# safety first
set -eu 

# commad-line args:
test $# -lt 2 && {
   echo "usage: $0 fastq_path contigs_path"
   echo "e.g. $0 /work/m4b_binning/assembly/data/HOW10/8777.1.112183.AATAGG.fastq.gz /work/m4b_binning/assembly/longer_contigs/contigs_longer_than_1500bp.fa"
   exit
}

fastq_path=$1
contigs_path=$2

threads=1
#fastq_path=/work/m4b_binning/assembly/data/HOW10/8777.1.112183.AAGCGA.fastq.gz
#fastq_path=/work/m4b_binning/assembly/data/HOW10/8777.1.112183.AATAGG.fastq.gz
#contigs_path=/work/m4b_binning/assembly/longer_contigs/contigs_longer_than_1500bp.fa

# parse 8777.1.112183.AAGCGA type string from path like ...data/HOW10/8777.1.112183.AAGCGA.fastq.gz
sample=$(echo $fastq_path | sed 's#.*/\(.*\).fastq.gz#\1#') 
echo "parsed $sample from $fastq_path"

# parse the contig_details and put these mappings in there
contigs_base=$(echo $contigs_path| sed 's#.*/\(.*\).fa.*#\1#')
echo "parsed $contigs_base from $contigs_path"

# make a custom folder name, so we can map to multiple references (e.g. different min contig length)
dir="map_to_$contigs_base"
mkdir -p $dir

# One folder per sample.
dir=$dir/$sample
mkdir -p $dir
cd $dir

# First check that the run isn't complete.  Should have .sorted.bam if so
fastq=$(basename "$fastq_path")
if [ -e $fastq.sorted.bam ]
then
	echo "already found $fastq.sorted.bam"
	echo "do not re-run this combination"
	exit 1
fi


# Want to keep a record of how long the script took.  Cat the time before & after to a file. 
t_file='timestamps.txt'
touch $t_file
echo "$(date)" >> $t_file 

# Move the reags fastq and the contigs fasta to the new folder
cp $fastq_path .
cp $contigs_path .
contigs=$(basename "$contigs_path")
echo "do stuff on $fastq"

# Index both files
if [ ! -e $sample.bwt ]
then 
	bwa index $contigs
	echo "done with bwa index"
	samtools faidx $contigs
	echo "done with bwa samtools faidx"
fi

## the first part maps the fastqgz to the bin fasta that pipes to samtools to create a binary that is fed to another samtools to sort by contig
cmd="bwa mem -t $threads $contigs $fastq | samtools view -h -b -S /dev/stdin | samtools sort -m 1000000000 -o $fastq.sorted.bam /dev/stdin"
echo $cmd
echo $cmd >> $t_file   # record the main bwa/samtools call 
eval $cmd

echo "$(date)" >> $t_file 

cd ..
