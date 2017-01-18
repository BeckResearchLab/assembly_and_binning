# sample fastq:
# /work/m4b_binning/assembly/data/HOW11/8777.2.112196.AGATAG.fastq.gz

RESULTS='sample_read_counts.tsv'
rm -f $RESULTS

for fq in `ls ../*/*.fastq.gz`
do
	fastq=`basename $fq`
	#echo $fq
	# put the filename in without a new line:
	echo -n "$fastq	" >> $RESULTS 
	
	# append the number of @ lines in the file:
	count_command="zcat $fq | grep '^@' | wc -l >> $RESULTS "
	echo "count command: $count_command"
	eval $count_command
	

	##count_command="zcat $fq | grep '^@' | wc -l | >> $RESULTS "
	#count_command="zcat $fq | grep '^@' | wc -l"
	#echo "count command: $count_command"
	##num_reads=`$count_command`
	##num_reads=$(`$count_command`)
	#echo "number of reads for $fastq: $num_reads"
	##eval `echo $fastq $num_reads >> $RESULTS` 
done

echo "done counting all reads in each fastq.gz file"
