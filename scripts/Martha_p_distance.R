




# Base R function


read_fasta_sequence <- function(fasta_file) {
  lines <- readLines(fasta_file)
  # Assuming the first line is the header and subsequent lines are the sequence
  sequence <- paste(lines[-1], collapse = "")
  return(sequence)
}




calc_p_distance <- function(path_to_ref, path_to_query) {
  # Read fasta file
  ref <- read_fasta_sequence(path_to_ref)
  query <- read_fasta_sequence(path_to_query)
  
  # Convert the sequences to character vectors
  ref_chars <- strsplit(ref, split = "")[[1]]
  query_chars <- strsplit(query, split = "")[[1]]
  
  # Count the differences (SNPs)
  snps <- sum(seq1_chars != seq2_chars)
  
  # Calculate p-distance
  p_dist <- snps/length(ref_chars)
  
  # output
  return(p_dist)
}

# Example usage

calc_p_distance(path_to_ref = "./data/tmp_ref_b99.fasta", path_to_query = "./data/tmp_query.fasta")




# end 








require(pacman)
pacman::p_load(Biostrings)

# fasta
ref <- readDNAStringSet("./data/tmp_ref_b99.fasta")
query <- readDNAStringSet("./data/tmp_query.fasta")

ref_header <- names(ref)

query_header <- names(query)
#####
read.fasta <- function(path_to_ref, path_to_seq){
  ref <- readDNAStringSet(path_to_ref)
  ref_header <- names(ref)
  query <- readDNAStringSet(path_to_seq)
  query_header <- names(query)
}


# SNPS (NO transitions vs tranversions)
count_nt_diff <- function(){
  
}



# SNPS-2 (distinguish transitions vs tranversions)




# P-distance
calc_p_distance <- function(){
  
}


# Jukes-Cantor (JC) model 




# Kimura 2-parameter (K2P) model (transitions vs tranversions)




# General Time Reversible (GTR) Model (complex!!!)





# Maximum Composite Likelihood method?


