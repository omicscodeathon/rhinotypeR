#' Visualize SNP differences against a reference (fast zoom & highlighting)
#'
#' @description
#' `SNPeek()` renders per–position nucleotide differences between an aligned
#' set of DNA sequences and a chosen reference using a cached pre-computation
#' for fast redraws. With no `xlim`/`center`+`window`, it plots the full
#' genome span and **all** sequences. Users can zoom regions and optionally
#' emphasize specific sequence IDs with(out) hiding the others.
#'
#' @param aln_set A [`Biostrings::DNAStringSet`] of aligned sequences
#'   (all sequences must have identical width).
#' @param ref_name Character; the sequence name in `aln_set` to use as
#'   reference. If `NULL`, the last sequence is used.
#' @param xlim Optional integer length-2 vector giving the genomic window,
#'   e.g. `c(100, 200)`. If omitted, the entire span is shown.
#' @param center,window Optionally specify a zoom by center position and
#'   window width (e.g., `center = 5012, window = 600`). Ignored if `xlim`
#'   is provided.
#' @param showLegend Logical; if `TRUE`, draws a legend for the nucleotide
#'   color map. Default: `FALSE`.
#' @param colorMapNT Named character vector mapping nucleotides to colors.
#'   Default: `c(A = "green3", T = "red2", C = "blue2", G = "gold")`.
#' @param colorFallback Color for any non-ATCG symbol (e.g., `N`, `-`).
#'   Default: `"grey50"`.
#' @param highlight_seqs Optional character vector of sequence IDs (row names)
#'   to emphasize (adds a translucent row band and thicker mismatch strokes).
#'   All sequences are still shown unless `show_only_highlighted = TRUE`.
#' @param highlight_bg Fill color for emphasized rows (use
#'   `grDevices::adjustcolor()` for transparency). Default:
#'   `adjustcolor("#FFE082", alpha.f = 0.5)`.
#' @param lwd_normal,lwd_highlight Line widths for mismatch segments in normal
#'   vs highlighted rows. Defaults: `1` and `2`.
#' @param seg_half_height Vertical half-height of each mismatch segment (row
#'   units). Default: `0.35`.
#' @param y_cex Expansion factor for y-axis (sequence labels). Default: `0.8`.
#' @param line_col Baseline axis/line color. Default: `"grey20"`.
#' @param show_only_highlighted Logical; if `TRUE`, display only rows whose
#'   IDs are in `highlight_seqs`. Default: `FALSE`.
#'
#' @details
#' Internally, `SNPeek()` builds a reusable cache of mismatch
#' positions and color indices versus the reference, then subsets to 
#' the requested window before plotting for speed.
#' The cache is returned invisibly so you can reuse it 
#' for interactive workflows without recomputation.
#'
#' @return Invisibly returns the precomputed cache (class `"SNPeekCache"`)
#'   used for plotting. The primary output is a base R plot.
#'
#' @seealso  `plotAA()` for the amino-acid analogue.
#'
#' @family visualization
#'
#' @examples
#' fasta_file <- system.file("extdata", "input_aln.fasta", package = "rhinotypeR")
#' aln <- Biostrings::readDNAStringSet(fasta_file)
#'
#' # Full span
#' SNPeek(aln)
#'
#' # Zoom + highlight
#' SNPeek(aln, xlim = c(100, 200),
#'        highlight_seqs = c("MT177780.1", "MT177798.1"),
#'        showLegend = TRUE)
#' 
#' @importFrom graphics par plot axis legend rect segments
#' @importFrom grDevices adjustcolor
#' @export
SNPeek <- function(aln_set,
                   ref_name = NULL,
                   xlim = NULL,
                   center = NULL, window = NULL,
                   showLegend = FALSE,
                   colorMapNT = c(A = "green3", T = "red2", C = "blue2", G = "gold"),
                   colorFallback = "grey50",
                   highlight_seqs = NULL,
                   highlight_bg = grDevices::adjustcolor("#FFE082", alpha.f = 0.5),
                   lwd_normal = 1,
                   lwd_highlight = 2,
                   seg_half_height = 0.35,
                   y_cex = 0.8,
                   line_col = "grey20",
                   show_only_highlighted = FALSE) {
  
  # Build cache (compute-once if you want to manage it outside SNPeek)
  cache <- compute_cache(
    aln_set = aln_set,
    ref_name = ref_name,
    colorMap = colorMapNT,
    colorFallback = colorFallback
  )
  
  # If user didn’t specify xlim/center+window → plot full genome span
  # (plot_window will default to full span when xlim is NULL)
  plot_window(cache = cache,
              xlim = xlim, center = center, window = window,
              showLegend = showLegend,
              highlight_seqs = highlight_seqs,
              highlight_bg = highlight_bg,
              lwd_normal = lwd_normal, lwd_highlight = lwd_highlight,
              seg_half_height = seg_half_height,
              y_cex = y_cex,
              line_col = line_col,
              show_only_highlighted = show_only_highlighted)
  
  invisible(cache) # return cache to reuse for rapid re-draws
}
