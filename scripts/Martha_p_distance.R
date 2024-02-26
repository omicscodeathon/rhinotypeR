#install.packages("devtools")
devtools::install_github("gtonkinhill/pairsnp-r")


library(pairsnp)
fasta.file.name <- system.file("extdata", "seqs.fa", package = "pairsnp")
sparse.data <- import_fasta_sparse(fasta.file.name)
d <- snp_dist(sparse.data)
s <- snp_similarity(sparse.data)


# start


require(pacman)
pacman::p_load(Biostrings)

# fasta
ref <- readDNAStringSet("./data/tmp_ref_b99.fasta")
query <- readDNAStringSet("./data/tmp_query.fasta")

ref_header <- names(ref)

query_header <- names(query)


read_fasta_sequence <- function(fasta_file) {
  lines <- readLines(fasta_file)
  # Assuming the first line is the header and subsequent lines are the sequence
  sequence <- paste(lines[-1], collapse = "")
  return(sequence)
}

# Example usage
sequence1 <- read_fasta_sequence("./data/tmp_ref_b99.fasta")
sequence2 <- read_fasta_sequence("./data/tmp_query.fasta")


# Convert the sequences to character vectors
seq1_chars <- strsplit(sequence1, split = "")[[1]]
seq2_chars <- strsplit(sequence2, split = "")[[1]]

# Count the differences (SNPs)
snps <- sum(seq1_chars != seq2_chars)

# Print the number of SNPs
print(snps)


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


