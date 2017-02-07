set -eu

fastq_files=`ls /work/m4b_binning/assembly/data/*OW*/*.fastq.gz` # `` same as $()
echo "fastq files found:"
# You don't alwys need to echo something in quotes.
echo "$fastq_files" | wc -w

# Doesn't work: it isn't splititng up the fastq files. 
#echo $fastq_files | parallel --jobs 20 map_reads.sh {} /work/m4b_binning/assembly/longer_contigs/contigs_longer_than_1500bp.fa

# test that worked:
#ls /work/m4b_binning/assembly/data/*OW*/*.fastq.gz | parallel --jobs 20 ./args.sh {} /work/m4b_binning/assembly/longer_contigs/contigs_longer_than_1500bp.fa

ls /work/m4b_binning/assembly/data/*OW*/*.fastq.gz | parallel --jobs 20 ./map_reads.sh {} /work/m4b_binning/assembly/longer_contigs/contigs_longer_than_1500bp.fa
#ls /work/m4b_binning/assembly/data/*OW*/*.fastq.gz | parallel --jobs 20 old/args.sh {} /work/m4b_binning/assembly/longer_contigs/contigs_longer_than_1500bp.fa
