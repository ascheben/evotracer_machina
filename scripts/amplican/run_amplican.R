library(amplican)

args <- commandArgs(trailingOnly = TRUE)
fastq_dir <- args[1]
config_filepath <- args[2]
output_dir <- args[3]

## for testing amplican function before integration into pipeline
#fastq_dir <- "sim_results_test/out_simulator_test/all_samples_fastq"
#config_filepath <- "sim_results_test/out_amplican_test/amplican_config_test.csv"
#output_dir <- "sim_results_test/out_amplican_test/"
 
#  run amplican
amplicanPipeline(config_filepath, 
                fastq_dir, 
                output_dir,
                knit_reports = FALSE,
                event_filter=TRUE)

# results of the analysis can be found at
#message(output_dir)

