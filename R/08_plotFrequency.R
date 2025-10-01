#' Plot Frequency of Assigned Rhinovirus Types
#'
#' Generates a barplot showing the frequency of assigned rhinovirus types.
#' Optionally displays a legend indicating species classification (A, B, C, Other).
#'
#' @param assignedTypesDF A data frame, typically the output of \code{\link{assignTypes}},
#'   containing at least the columns \code{query} and \code{assignedType}.
#' @param showLegend Logical; if \code{TRUE}, adds a legend mapping species
#'   (A, B, C, Other) to colors. Default is \code{FALSE}.
#' @param sort Logical; if \code{TRUE}, bars are sorted descending by count. Default \code{FALSE}.
#' @param las Integer; axis label style for type names (defaults to 2 = perpendicular for readability).
#' @param cex_names Numeric; character expansion for x-axis names. Default 0.8.
#'
#' @details
#' Types are grouped into species via the first letter \emph{of the type token}
#' after stripping common prefixes like \code{"RV-"} or \code{"RV"}.
#' Any \code{"unassigned"} (any case) is labeled as species \code{"Other"}.
#'
#' Colors:
#' \itemize{
#'   \item A = blue
#'   \item B = red
#'   \item C = green
#'   \item Other = grey
#' }
#'
#' @return Invisibly returns a data frame with columns:
#'   \code{assignedType}, \code{count}, \code{species}.
#'
#' @family visualization
#' @seealso \code{\link{assignTypes}}
#'
#' @examples
#' test <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")
#' fastaD <- Biostrings::readDNAStringSet(test)
#' assigned <- try(assignTypes(fastaD), silent = TRUE)
#' if (!inherits(assigned, "try-error")) {
#'   plotFrequency(assigned, showLegend = TRUE)
#' }
#'
#' @importFrom graphics barplot legend
#' @export
plotFrequency <- function(assignedTypesDF,
                          showLegend = FALSE,
                          sort = FALSE,
                          las = 2,
                          cex_names = 0.8) {
  # basic checks
  if (!is.data.frame(assignedTypesDF)) {
    stop("assignedTypesDF must be a data.frame (output of assignTypes).")
  }
  req <- c("query", "assignedType")
  missing_cols <- setdiff(req, names(assignedTypesDF))
  if (length(missing_cols)) {
    stop("assignedTypesDF is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # normalize assignedType
  atype <- as.character(assignedTypesDF$assignedType)
  if (length(atype) == 0L) stop("No rows in assignedTypesDF.")
  
  # count per assigned type
  types_counts <- as.data.frame(table(assignedType = atype), stringsAsFactors = FALSE)
  names(types_counts) <- c("assignedType", "count")
  
  # species mapping helper:
  # strip common prefixes (RV-, RV_, RV), then take first letter; unassigned -> Other
  norm_first_letter <- function(x) {
    ifelse(grepl("^\\s*unassigned\\s*$", x, ignore.case = TRUE),
           "Other",
           {
             # remove RV prefix variants and delimiters
             x2 <- gsub("^(RV[-_ ]?|RV)", "", x, ignore.case = TRUE)
             fl <- toupper(substr(x2, 1, 1))
             ifelse(fl %in% c("A","B","C"), fl, "Other")
           })
  }
  
  types_counts$species <- vapply(types_counts$assignedType, norm_first_letter, character(1))
  
  # colors
  color_map <- c(A = "blue", B = "red", C = "green", Other = "grey")
  bar_cols <- unname(color_map[types_counts$species])
  
  # optional sorting
  if (isTRUE(sort)) {
    o <- order(types_counts$count, decreasing = TRUE)
    types_counts <- types_counts[o, , drop = FALSE]
    bar_cols <- bar_cols[o]
  }
  
  # plot
  bp <- barplot(
    height = types_counts$count,
    names.arg = types_counts$assignedType,
    col = bar_cols,
    main = "Frequency of Assigned Types",
    xlab = "RV Type",
    ylab = "Count",
    las = las,
    cex.names = cex_names
  )
  
  if (isTRUE(showLegend)) {
    legend("topright",
           legend = names(color_map),
           fill   = unname(color_map),
           title  = "Species",
           bty    = "n",
           inset  = 0.02,
           horiz  = FALSE)
  }
  
  invisible(types_counts)
}
