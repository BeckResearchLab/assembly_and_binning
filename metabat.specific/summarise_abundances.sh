source activate py3
set -eu

# parse the bins to make a tsv containing contig <--> bin mappings
python ../abundances/bin_contig_analysis.py -input . -output bin_contig_mappings.tsv

# aggregate abundances for each contig, roll up by bin.  
python ../abundances/aggregate.py -contig_path ./bin_contig_mappings.tsv -depth_path ./depth.txt -sample_info_dir ../data/sample_info/
