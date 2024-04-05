# Download prototype strains

# @Ruth, the data will now be a system.file, package = rhinotypeR, 
#so adjust the 'from' paths accordingly

getPrototypeSeqs <- function(destinationFolder){
  
  # copy files
  file.copy(from = "./data/RVAPrototypeAligned.fasta", 
            to = file.path(destinationFolder,"RVARefs.fasta"), overwrite = TRUE)
  file.copy(from = "./data/RVBPrototypeAligned.fasta", 
            to = file.path(destinationFolder,"RVBRefs.fasta"), overwrite = TRUE)
  file.copy(from = "./data/RVCPrototypeAligned.fasta", 
            to = file.path(destinationFolder,"RVCRefs.fasta"), overwrite = TRUE)
  
  print(paste0("The reference sequences have been downloaded to ", destinationFolder))
  
}

# Example usage (download into the output folder)
getPrototypeSeqs(destinationFolder = "./output")

