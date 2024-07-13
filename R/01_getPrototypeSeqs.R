# Download prototype strains

getPrototypeSeqs <- function(destinationFolder){
  
  ref <- system.file("extdata", "prototypes.fasta", package = "rhinotypeR")
  # copy files
  file.copy(from = ref,
            to = file.path(destinationFolder,"RVRefs.fasta"), overwrite = TRUE)
  
  message("The reference sequences have been downloaded to ", destinationFolder)
  
}
