bam_files=`ls ../*.bam`
# echo $bam_files

mkdir -p idxstat_results

files_done=0
for bam in $bam_files
do
	bam_path=`realpath $bam`
	echo "loop $(($files_done + 1)): index and get idxstats for $bam_path"
	# index the .bam file 
	samtools index $bam

	# do idxstats.  Seems to want the .bai in the same dir. 
	idxstats_out_path=`basename $bam`".idxstats"
	idxcommand="samtools idxstats $bam > ./idxstat_results/$idxstats_out_path"
	echo "idxcommand: $idxcommand"
	eval $idxcommand

	# remove .bai file (they are quick to make, so no biggie)
	#rm $bam.bai

	# add count for completed call 
	files_done=$((files_done + 1))
done
