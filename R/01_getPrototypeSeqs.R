#' Download rhinovirus prototype strains (optional)
#'
#' Copies the packaged rhinovirus prototype strains
#' (\file{prototypes.fasta}) from \pkg{rhinotypeR} into a user-specified
#' directory as \file{RVRefs.fasta}. Downloading to local storage is **not
#' required** to use the package.
#'
#' Users have two equivalent workflows:
#' \enumerate{
#'   \item \strong{In-R workflow (no download):} use
#'         \code{\link{alignToRefs}} to align your sequences against the
#'         packaged reference prototypes directly in R, then run
#'         \code{\link{assignTypes}} on the resulting alignment.
#'   \item \strong{External-tools workflow (with download):} use
#'         \code{getPrototypeSeqs} to save \file{RVRefs.fasta} locally,
#'         \emph{combine it with your new sequences}, and perform alignment
#'         using your preferred external tool. You can then bring the aligned
#'         FASTA back into R for genotype assignment.
#' }
#' 
#' @param destinationFolder \code{character(1)}. Path to an existing directory
#'   where the prototype FASTA will be written. The folder must already exist.
#' @param overwrite \code{logical(1)}. Whether to overwrite an existing
#'   \file{RVRefs.fasta} in the destination directory. Defaults to \code{TRUE}.
#'
#' @details
#' Internally, this function uses \code{\link[base]{system.file}} to locate the
#' packaged prototype file and \code{\link[base]{file.copy}} to copy it into the
#' specified directory. If you are following the in-R workflow, you can skip this
#' function entirely and proceed with \code{\link{alignToRefs}} followed by
#' \code{\link{assignTypes}}.
#'
#' @return
#'  Prints a message noting the destination directory.
#'
#' @author
#' Martha Luka, Ruth Nanjala, Wafaa Rashed, Winfred Gatua, Olaitan Awe
#'
#' @seealso
#' \code{\link{assignTypes}}, \code{\link{alignToRefs}}
#'
#' @examples
#' if (interactive()) {
#'   # --- In-R workflow (no download) ---
#'   # aln <- alignToRefs(query_fasta = "my_new_sequences.fasta", method = "Muscle")
#'   # typesDF <- assignTypes(aln)
#'
#'   # --- External-tools workflow (with download) ---
#'   dest_dir <- tempdir() # specify a destination directory
#'   getPrototypeSeqs(destinationFolder = dest_dir)  # writes RVRefs.fasta
#'   # Now combine RVRefs.fasta with your new sequences and align using
#'   # external tools (e.g., MAFFT, MUSCLE). Then read the aligned FASTA
#'   # back into R and proceed with assignTypes().
#'   list.files(dest_dir)
#' }
#'
#' @export
getPrototypeSeqs <- function(destinationFolder, overwrite = TRUE) {
  # Check if folder exists
  if (!dir.exists(destinationFolder)) 
    stop("destinationFolder does not exist: ", destinationFolder)
  
  # Gets the prototypes
  ref <- system.file("extdata", "prototypes.fasta", package = "rhinotypeR")
  
  # Copy file to destination
  file.copy(from = ref,
            to = file.path(destinationFolder, "RVRefs.fasta"),
            overwrite = overwrite)
  
  message("The reference sequences have been downloaded to ", destinationFolder)
}
