source activate py2  # BioPython required

python shorten_contig_names.py



# [ec2-user@ip-10-0-0-158 genomes]$ head -n 1 nucleotide/Methylophilus_sp_Q8.fa
# >GQ52DRAFT_scf7180000000002_quiver_dupTrim_7654.1       >2585443431 GQ52DRAFT_scf7180000000002_quiver_dupTrim_7654.1


# Dumps tons of stuff to the console for 5 min!
#for genome in ./nucleotide/*.fa
#do 
#	echo "sed -re 's/(_length)[^=]*$/\1/' ${genome}" 
#	sed -re 's/(_length)[^=]*$/\1/' ${genome}
#done
