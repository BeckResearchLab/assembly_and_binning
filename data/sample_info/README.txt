sample_info.tsv which had mappings from sample names to week, o2, replicate were found to be missing 5 rows
- this is true, also, for the copies on waffle.  I renamed them to signify they are broken.
- Found a verified copy on J's laptop:  /Users/janet/Dropbox/meta4/reference_data/sample_meta_info.tsv 
	- double checked all the fields (and the patterns between them) by hand.

See /work/jmatsen/170110_parse_meta4_sample_names

(py3)waffle:170110_parse_meta4_sample_names jmatsen$ ls
gather_sample_name_mappings.py  meta4_sample_names--cryptic_to_sample_number.tsv  meta4_sample_names--simple_only.tsv  meta4_sample_names.tsv
(py3)waffle:170110_parse_meta4_sample_names jmatsen$ scp -i ~/.ssh/janet_matsen.pem meta4_sample_names--cryptic_to_sample_number.tsv ec2-user@35.165.146.147:/work/m4b_binning/assembly/data/sample_info
meta4_sample_names--cryptic_to_sample_number.tsv                                                                                                                                                                                100% 6324     6.2KB/s   00:00
(py3)waffle:170110_parse_meta4_sample_names jmatsen$ ls
gather_sample_name_mappings.py  meta4_sample_names--cryptic_to_sample_number.tsv  meta4_sample_names--simple_only.tsv  meta4_sample_names.tsv
(py3)waffle:170110_parse_meta4_sample_names jmatsen$ pwd
/work/jmatsen/170110_parse_meta4_sample_names

