waffle:individual_isolates (master) jmatsen$ ls fastas/
Ancylobacter_sp_FA202.fa  Flavobacterium_sp_83.fa       Methylobacterium_sp_77.fa             Methylomonas_sp_MK1.fa            Methylopila_sp_M107.fa         Methylotenera_sp_1P_1.fa         Methylotenera_versatilis_7.fa            Mycobacterium_sp_155.fa
Arthrobacter_sp_31Y.fa    Flavobacterium_sp_Fl.fa       Methylobacterium_sp_88A.fa            Methylophilaceae_bacterium_11.fa  Methylosarcina_lacus_LW14.fa   Methylotenera_sp_73s.fa          Methyloversatilis_sp_FAM1.fa             Paracoccus_sp_N5.fa
Arthrobacter_sp_35W.fa    Hoeflea_sp_108.fa             Methylobacter_tundripaludum_21_22.fa  Methylophilus_sp_1.fa             Methylosinus_sp_LW3.fa         Methylotenera_sp_G11.fa          Methyloversatilis_sp_RZ18-153.fa         Pseudomonas_sp_11_12A.fa
Arthrobacter_sp_MA-N2.fa  Hyphomicrobium_sp_802.fa      Methylobacter_tundripaludum_31_32.fa  Methylophilus_sp_42.fa            Methylosinus_sp_LW4.fa         Methylotenera_sp_L2L1.fa         Methyloversatilis_universalis_Fam500.fa  Xanthobacteraceae_bacterium_501b.fa
Bacillus_sp_37MA.fa       Hyphomicrobium_sp_99.fa       Methylocystis_sp_LW5.fa               Methylophilus_sp_5.fa             Methylosinus_sp_PW1.fa         Methylotenera_sp_N17.fa          Methyloversatilis_universalis_FAM5.fa    Xanthobacter_sp_126.fa
Bacillus_sp_72.fa         Janthinobacterium_sp_RA13.fa  Methylomonas_sp_11b.fa                Methylophilus_sp_Q8.fa            Methylotenera_mobilis_13.fa    Methylotenera_versatilis_301.fa  Methylovorus_glucosetrophus_SIP3-4.fa    Xanthobacter_sp_91.fa
Bosea_sp_117.fa           Methylobacterium_sp_10.fa     Methylomonas_sp_LW13.fa               Methylopila_sp_73B.fa             Methylotenera_mobilis_JLW8.fa  Methylotenera_versatilis_79.fa   Mycobacterium_sp_141.fa


waffle:individual_isolates (master) jmatsen$ pwd
/dacb/meta4_iso/data/individual_isolates
waffle:individual_isolates (master) jmatsen$ scp -i ~/.ssh/janet_matsen.pem ./fastas/*.fa ec2-user@35.165.146.147:/work/m4b_binning/assembly/isolate_Fauzi_stats/isolate/genomes/nucleotide
