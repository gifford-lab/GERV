mcc -mv kmersnp_module -a preprocess/kmer_mat/correct_for_resol.m preprocess/baseline/parseBaseline.m preprocess/seq_data/findValidationCutoff.m snp-kmer-scoring/*.m postprocess/combineKSS/combine_snp.m
#cd snp-kmer-scoring
#mcc -mv main_baseline.m -a generateKmer.m aggregateKmer.m exp_errorL2_BL.m
#cd ../
