library(Biostrings)
library(tidyverse)

calculate_alignment_error <- function(target_seq, pool_seqs) {
  target_length <- nchar(target_seq)
  alignment_errors <- c()

  for (i in 1:length(pool_seqs)) {
    pool_seq <- as.character(pool_seqs[i])
    
    # Perform Needleman-Wunsch alignment
    sigma <- nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = TRUE)
    alignment <- pairwiseAlignment(DNAString(target_seq), DNAString(pool_seq), substitutionMatrix = sigma, gapOpening = -1, gapExtension = 0)
    
    # Calculate alignment error as the number of gaps and mismatches
    alignment_error <- target_length - alignment@score
    alignment_errors <- c(alignment_errors, alignment_error)
  }

  best_match_error <- min(abs(alignment_errors))
  return(list(best_match_error = best_match_error, all_errors = alignment_errors))
}


args <- commandArgs(trailingOnly = TRUE)
pool_sequences_file <- args[1]
df_file <- args[2]
out_dir <- args[3]

# Read the fasta file and extract the pool sequences
pool_sequences <- readDNAStringSet(pool_sequences_file)

# Read the CSV file and extract the target sequences
df <- read_csv(df_file)
df <- df %>% distinct(asv_names, .keep_all = TRUE)
target_indices <- df$asv_names
target_sequences <- df$seq

no_alignment_sequences <- NULL

# Perform sequence alignment and calculate alignment error
alignment_errors <- map2(target_sequences, target_indices, ~{
  target_seq <- as.character(.x)
  errors <- calculate_alignment_error(target_seq, pool_sequences)

  # Get the pool sequences that do not have an alignment error equal to the best error
  pool_sequences_no_alignment <- pool_sequences[which(errors$all_errors != errors$best_match_error)]
  
  # Initialize no_alignment_sequences with all pool sequences on the first occurrence
  if (is.null(no_alignment_sequences)) {
    no_alignment_sequences <<- pool_sequences_no_alignment
  } else {
    # Retain only the sequences found in both lists
    no_alignment_sequences <<- intersect(no_alignment_sequences, as.character(pool_sequences_no_alignment))
  }

  data.frame(index = .y, 
              best_alignment_error = errors$best_match_error,
              count_best_alignment_error = sum(errors$all_errors == errors$best_match_error),
              all_alignment_errors = paste(errors$all_errors, collapse = " "), 
              stringsAsFactors = FALSE)
})

# Create a dataframe with the alignment errors
alignment_df <- do.call(rbind, alignment_errors)

# Save the alignment errors to a CSV file
output_file <- file.path(out_dir, "alignment_errors.csv")
write_csv(alignment_df, output_file)

# Save the pool barcode sequences without a best alignment to a txt file
no_alignment_sequences_file <- file.path(out_dir, "barcode_no_ASV_alignment_sequences.txt")
writeLines(as.character(no_alignment_sequences), no_alignment_sequences_file)

