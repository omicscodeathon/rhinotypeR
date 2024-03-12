

# Base R functions


# Function 1 
# Read fasta
read_fasta_sequence <- function(fasta_file) {
  lines <- readLines(fasta_file)
  sequence <- paste(lines[-1], collapse = "")
  return(toupper(sequence)) # Convert sequence to uppercase
}

# Example usage
read_fasta_sequence("./data/tmp_query.fasta")
read_fasta_sequence("./data/tmp_ref_b99.fasta")

# Function 2
# Compare sequence lengths
compare_lengths <- function(path_to_ref, path_to_query){
  ref <- read_fasta_sequence(path_to_ref)
  query <- read_fasta_sequence(path_to_query)
  
  # Convert the sequences to character vectors
  ref_length <- length(unlist(strsplit(ref, "")))
  query_length <- length(unlist(strsplit(query, "")))
  
  # Return both lengths for further use
  return(list(ref_length = ref_length, query_length = query_length))
}
# Example usage
compare_lengths(path_to_ref = "./data/tmp_ref_b99.fasta", path_to_query="./data/tmp_query.fasta")

# Function 3
# Function to handle gaps or non-ATCG characters in sequences
handle_gaps <- function(ref_seq, query_seq, option = "delete") {
  # Convert sequences to character vectors
  ref_chars <- strsplit(toupper(ref_seq), "")[[1]]
  query_chars <- strsplit(toupper(query_seq), "")[[1]]
  
  # Identify positions with non-ATCG characters
  valid_positions <- which((ref_chars %in% c("A", "T", "C", "G")) & (query_chars %in% c("A", "T", "C", "G")))
  
  if (option == "delete") {
    # Keep only positions with ATCG characters in both sequences
    ref_chars <- ref_chars[valid_positions]
    query_chars <- query_chars[valid_positions]
  } # If option is "ignore", we do nothing to the sequences
  
  # Reassemble sequences
  ref_seq_cleaned <- paste(ref_chars, collapse = "")
  query_seq_cleaned <- paste(query_chars, collapse = "")
  
  return(list(ref_seq = ref_seq_cleaned, query_seq = query_seq_cleaned))
}



# Function 4
# Count SNPs
count_SNPs <- function(path_to_ref, path_to_query){
  
  # Read fasta file
  ref <- read_fasta_sequence(path_to_ref)
  query <- read_fasta_sequence(path_to_query)
  
  # Convert the sequences to character vectors
  ref_chars <- strsplit(ref, split = "")[[1]]
  query_chars <- strsplit(query, split = "")[[1]]
  
  # Count the differences (SNPs)
  snps <- sum(ref_chars != query_chars)
  
  # output
  return(snps)
}
# Example usage
count_SNPs(path_to_ref = "./data/tmp_ref_b99.fasta", path_to_query="./data/tmp_query.fasta")




# Function 5
# p-distance
calc_p_distance <- function(path_to_ref, path_to_query) {
  
  # Count the differences (SNPs)
  snps <- count_SNPs(path_to_ref, path_to_query)
  
  # lengths
  lengths <- compare_lengths(path_to_ref, path_to_query)
  ref_length <- lengths$ref_length
  
  # Calculate p-distance
  p_dist <- snps/ref_length
  
  # output
  return(p_dist)
}

# Example usage
calc_p_distance(path_to_ref = "./data/tmp_ref_b99.fasta", path_to_query = "./data/tmp_query.fasta")


# Jukes Cantor

# Function to calculate Jukes-Cantor genetic distance
calc_jukes_cantor_distance <- function(path_to_ref, path_to_query) {
  # Calculate p-distance
  p_dist <- calc_p_distance(path_to_ref, path_to_query)
  
  # Calculate Jukes-Cantor genetic distance
  jc_dist <- -3/4 * log(1 - 4/3 * p_dist)
  
  # Return the Jukes-Cantor genetic distance
  return(jc_dist)
}

# Example usage
calc_jukes_cantor_distance(path_to_ref = "./data/tmp_ref_b99.fasta", path_to_query = "./data/tmp_query.fasta")


# Kimura 2 parameter (transitions vs transversions)

## Function to calculate proportions of transitions and transversions
calculate_transitions_transversions <- function(ref_chars, query_chars) {
  transitions <- sum((ref_chars == 'A' & query_chars == 'G') |
                       (ref_chars == 'G' & query_chars == 'A') |
                       (ref_chars == 'C' & query_chars == 'T') |
                       (ref_chars == 'T' & query_chars == 'C'))
  
  transversions <- sum((ref_chars == 'A' & query_chars == 'C') |
                         (ref_chars == 'A' & query_chars == 'T') |
                         (ref_chars == 'G' & query_chars == 'C') |
                         (ref_chars == 'G' & query_chars == 'T') |
                         (ref_chars == 'C' & query_chars == 'A') |
                         (ref_chars == 'C' & query_chars == 'G') |
                         (ref_chars == 'T' & query_chars == 'A') |
                         (ref_chars == 'T' & query_chars == 'G'))
  
  total_sites <- length(ref_chars)
  P <- transitions / total_sites
  Q <- transversions / total_sites
  
  return(list(P = P, Q = Q))
}

## Function to calculate Kimura 2-parameter genetic distance
calc_kimura_2p_distance <- function(path_to_ref, path_to_query) {
  
  # Read fasta files and ensure sequences are in uppercase
  ref <- read_fasta_sequence(path_to_ref)
  query <- read_fasta_sequence(path_to_query)
  
  # Convert the sequences to character vectors
  ref_chars <- strsplit(ref, split = "")[[1]]
  query_chars <- strsplit(query, split = "")[[1]]
  
  # Calculate proportions of transitions and transversions
  tt <- calculate_transitions_transversions(ref_chars, query_chars)
  P <- tt$P
  Q <- tt$Q
  
  # Calculate Kimura 2-parameter genetic distance
  K2P_distance <- -0.5 * log((1 - 2*P - Q) * sqrt(1 - 2*Q))
  
  # Return the Kimura 2-parameter genetic distance
  return(K2P_distance)
}

# Example usage
calc_kimura_2p_distance(path_to_ref = "./data/tmp_ref_b99.fasta", path_to_query = "./data/tmp_query.fasta")


# calculate Tamura 3-parameter genetic distance
calc_tamura_3p_distance <- function(path_to_ref, path_to_query) {
  # read_fasta_sequence 
  ref <- toupper(read_fasta_sequence(path_to_ref))
  query <- toupper(read_fasta_sequence(path_to_query))
  
  # Convert sequences to character vectors
  ref_chars <- strsplit(ref, "")[[1]]
  query_chars <- strsplit(query, "")[[1]]
  
  # Calculate base frequencies for both sequences
  gc_content_ref <- sum(ref_chars %in% c("G", "C")) / length(ref_chars)
  gc_content_query <- sum(query_chars %in% c("G", "C")) / length(query_chars)
  
  # Calculate transitions and transversions
  transitions <- sum((ref_chars == "A" & query_chars == "G") | (ref_chars == "G" & query_chars == "A") |
                       (ref_chars == "C" & query_chars == "T") | (ref_chars == "T" & query_chars == "C"))
  transversions <- sum((ref_chars != query_chars) & !((ref_chars == "A" & query_chars == "G") | (ref_chars == "G" & query_chars == "A") |
                                                        (ref_chars == "C" & query_chars == "T") | (ref_chars == "T" & query_chars == "C")))
  
  positions_scored <- sum(ref_chars != "-" & query_chars != "-" & !is.na(ref_chars) & !is.na(query_chars))
  P <- transitions / positions_scored
  Q <- transversions / positions_scored
  
  theta1 <- gc_content_ref
  theta2 <- gc_content_query
  C <- theta1 + theta2 - 2 * theta1 * theta2
  
  # Tamura distance calculation
  if ((1 - P/C - Q) > 0 && (1 - 2*Q) > 0) {
    distance <- -C * log(1 - P/C - Q) - 0.5 * (1 - C) * log(1 - 2*Q)
  } else {
    distance <- NA  # Distance cannot be calculated under these conditions
  }
  
  return(distance)
}


# Example usage
calc_tamura_3p_distance(path_to_ref = "./data/tmp_ref_b99.fasta", path_to_query = "./data/tmp_query.fasta")


# Tamura-Nei model  (complex!!!)

# General Time Reversible (GTR) Model (complex!!!)

# Maximum Composite Likelihood method? (possibly too complex for regular scripting!!!)


