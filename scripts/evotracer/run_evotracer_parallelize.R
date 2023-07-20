# development branch
###### libraries load ######
library("EvoTraceR")
###### Running from Copernicus ######

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[2]
trimmomatic_path <- args[3]
flash_path <- args[4]
devtools::load_all(args[5])

#devtools::load_all("/local/storage/no-backup/staklins-scratch/EvoTraceR-parallelize")
#input_dir <- "/local/storage/no-backup/staklins-scratch/evotracer_machina/fastq"
#output_dir <- "/local/storage/no-backup/staklins-scratch/evotracer_machina/out_parallelize"
#trimmomatic_path <- "/home/staklins/miniconda3/envs/r_env/share/trimmomatic-0.39-2/trimmomatic.jar"
#flash_path <- "/home/staklins/miniconda3/envs/r_env/bin/flash"

EvoTraceR_object <- 
  initialize_EvoTraceR(
    input_dir = input_dir,
    output_dir = output_dir,
    map_file_sample = NULL,
    trimmomatic_path = trimmomatic_path,
    flash_path = flash_path)

start_time <- Sys.time()
###### filtering and calculating asv statistics ######
EvoTraceR_object <- 
  asv_analysis(EvoTraceR_object = EvoTraceR_object,
               ### Barcode Type ###
               ### BC10v0 ### 
               ref_name = "BC10v0", 
               ref_seq = "TCTACACGCGCGTTCAACCGAGGAAAACTACACACACGTTCAACCACGGTTTTTTACACACGCATTCAACCACGGACTGCTACACACGCACTCAACCGTGGATATTTACATACTCGTTCAACCGTGGATTGTTACACCCGCGTTCAACCAGGGTCAGATACACCCACGTTCAACCGTGGTACTATACTCGGGCATTCAACCGCGGCTTTCTGCACACGCCTACAACCGCGGAACTATACACGTGCATTCACCCGTGGATC", 
               ref_flank_left = "^TCTAC", # starting sequence of the barcode based on "ref_seq"  
               ref_flank_right = "CCCGTGGATC$", # ending sequence of the barcode based on "ref_seq"  
               ### general barcode features ###
               flanking_filtering = "either", # use both flanking sequences: c("left", "right", "both")
               ref_cut_sites = c(17, 43, 69, 95, 121, 147, 173, 199, 225, 251), # cut sites for Cas9 based on "ref_seq" 
               ref_border_sites = c(1, 26, 52, 78, 104, 130, 156, 182, 208, 234), # border for marking site (guide=20, PAM=3, spacer=3) = 26 x bp
               output_figures = TRUE,
               asv_count_cutoff = 3, # minimum number of ASVs to be counted; decided on: 03/25/22
               # pair-wise alignment parameters
               pwa_type = "global", # based on AmpliCan (global = Needleman-Wunsch)
               pwa_gapOpening = -25, # based on AmpliCan: -25
               pwa_gapExtension = 0, # based on AmpliCan: 0
               pwa_match = 15, # based on AmpliCan: 15
               pwa_mismatch = -4, # based on AmpliCan: -4
               cleaning_window = c(3,3), # cleaning window +/-3 from Cas9 cutting size 
               batch_size = 100,
               cores = parallel::detectCores()) 
end_time <- Sys.time()
end_time - start_time

prefilt <- EvoTraceR_object$asv_prefilter

seqtab <- EvoTraceR_object$clean_asv_dataframe

#?asv_analysis()
#?infer_phylogeny()
###### marking indexing and stat ######
EvoTraceR_object <- 
  analyse_mutations(EvoTraceR_object = EvoTraceR_object) # -> marking_analysis()

###### infer phylogeny ######
# Can this be multithreaded?
# phylogeny can be built based on: "del_ins", "smooth_del_ins", "del", "smooth_del"
EvoTraceR_object <- 
  infer_phylogeny(EvoTraceR_object = EvoTraceR_object, mutations_use = "del_ins")

###### plot treestat_bc ######
EvoTraceR_object <- 
  create_df_summary(EvoTraceR_object) # tree_asv_freq_plot

##### plot separated treestat_bc ######
# npt in the package
# plot_cluster_summary(EvoTraceR_object)

###### SAVE ALL DATA ######
#save.image(paste0(output_dir, "group_id_select", ".RData"))
###### SAVE ALL DATA ######
#save.image("MMUS1469_2023-03-23_Armin.RData")
