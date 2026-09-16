
# args for running interactively in R
args <- list()
args$alg_file <- '../sup_tables/TableS4_Merian_element_definitions.tsv'
args$busco1 <- '../busco_results/fol_ang/run_arthropoda_odb10/full_table.tsv'
args$busco2 <-'../busco_results/smi_aqu/run_arthropoda_odb10/full_table.tsv'
#args$busco3 <- 'test_data/Inachis_io.tsv'
args$chrom1 <- '../busco_results/fol_ang/fol_ang_info.tsv'
args$chrom2 <- '../busco_results/smi_aqu/smi_aqu_info.tsv'
#args$chrom3 <- 'test_data/Inachis_io_info.tsv'
args$alpha <- 0
args$output_prefix <- 'smi_aqu_fol_ang'
args$busco_list <- c('../busco_results/smi_aqu/run_arthropoda_odb10/full_table.tsv', '../busco_results/fol_ang/run_arthropoda_odb10/full_table.tsv')
args$chrom_list <- c('../busco_results/smi_aqu/smi_aqu_info.tsv', '../busco_results/fol_ang/fol_ang_info.tsv')
