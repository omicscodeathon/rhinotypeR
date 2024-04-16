# Download prototype strains

getPrototypeSeqs <- function(destinationFolder){
  
  # copy files
  file.copy(from = RVAPrototype,
            to = file.path(destinationFolder,"RVARefs.fasta"), overwrite = TRUE)
  file.copy(from = RVBPrototype,
            to = file.path(destinationFolder,"RVBRefs.fasta"), overwrite = TRUE)
  file.copy(from = RVCPrototype,
            to = file.path(destinationFolder,"RVCRefs.fasta"), overwrite = TRUE)
  
  print(paste0("The reference sequences have been downloaded to ", destinationFolder))
  
}

# Example usage (download into the output folder)
getPrototypeSeqs(destinationFolder = "./output")

