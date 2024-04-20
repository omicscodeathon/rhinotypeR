# Download prototype strains

getPrototypeSeqs <- function(destinationFolder){
  
  # copy files
  file.copy(from = "./data/prototypes_aln.fasta",
            to = file.path(destinationFolder,"RVRefs.fasta"), overwrite = TRUE)

  print(paste0("The reference sequences have been downloaded to ", destinationFolder))
  
}

# Example usage (download into the output folder)
getPrototypeSeqs(destinationFolder = "./output")

