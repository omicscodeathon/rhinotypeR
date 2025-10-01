#' Visualize amino acid differences against a reference (fast zoom & highlighting)
#'
#' @description
#' `plotAA()` plots per–position amino acid differences between an aligned
#' set of protein sequences and a chosen reference. It reuses the same cache and
#' fast redraw machinery as [SNPeek()], with an amino acid color map grouping
#' residues into biochemical classes:
#'
#' * Positively charged (Arg, His, Lys) – red
#' * Negatively charged (Asp, Glu) – blue
#' * Polar uncharged (Ser, Thr, Asn, Gln) – green
#' * Nonpolar / aromatic / special (all others) – gold
#' * Other/unknown (X, gaps) – grey
#'
#' By default (no `xlim`/`center`+`window`), the function plots the full
#' alignment span and all sequences. Users can zoom to regions of interest and
#' optionally emphasize specific sequence IDs without hiding the others.
#'
#' @param aln_aa_set A `Biostrings::AAStringSet` of aligned amino acid
#'   sequences (all sequences must be the same width).
#' @param ref_name Character scalar; the sequence name in `aln_aa_set` to use as
#'   reference. If `NULL`, the last sequence is used.
#' @param xlim Optional integer length-2 vector giving the protein window to
#'   display, e.g. `c(5, 20)`. If omitted, the entire span is shown.
#' @param center,window Optional integers; alternatively specify a zoom by
#'   center position and window width (e.g., `center = 150, window = 50`).
#'   Ignored if `xlim` is provided.
#' @param showLegend Logical; if `TRUE`, draws a compact legend for the five
#'   amino acid classes. Default: `FALSE`.
#' @param colorMapAA Named character vector mapping amino acids to colors. The
#'   default assigns one color per residue class (see Details).
#' @param colorFallback Color used for any symbol not present in `colorMapAA`
#'   (e.g., `X`, `-`). Default: `"grey50"`.
#' @param highlight_seqs Optional character vector of sequence IDs (row names)
#'   to visually emphasize in the plot (adds a translucent row band and thicker
#'   mismatch strokes). Sequences are still shown unless
#'   `show_only_highlighted = TRUE`.
#' @param highlight_bg Fill color for emphasized rows (semi-transparent
#'   recommended). Default: translucent amber.
#' @param lwd_normal,lwd_highlight Line widths for mismatch segments in normal
#'   vs highlighted rows. Defaults: `1` and `2`.
#' @param seg_half_height Vertical half-height of each mismatch segment (in row
#'   units). Default: `0.35`.
#' @param y_cex Expansion factor for y-axis (sequence labels). Default: `0.8`.
#' @param line_col Baseline axis/line color. Default: `"grey20"`.
#' @param show_only_highlighted Logical; if `TRUE`, display only the rows whose
#'   IDs are listed in `highlight_seqs`. Default: `FALSE`.
#'
#' @details
#' Internally, `plotAA()` calls [compute_cache()] with an amino acid–specific
#' color map, then passes the result to [plot_window()] for fast drawing. The
#' returned cache can be reused for repeated zooms without recomputation.
#'
#' When `showLegend = TRUE`, the legend shows only five biochemical classes
#' (positively charged, negatively charged, polar uncharged, nonpolar/aromatic,
#' and other).
#'
#' @return Invisibly returns the precomputed cache (class `"SNPeekCache"`)
#'   used for plotting. The primary output is a base R plot drawn to the
#'   active device.
#'
#' @seealso [SNPeek()] for the nucleotide version; [compute_cache()],
#'   [plot_window()] for internals.
#'
#' @family visualization
#'
#' @examples
#' # Create a protein alignment 
#' aa <- Biostrings::AAStringSet(c(
#'   Ref = "MTEYKLVVVGYKL",
#'   S1  = "MTEYKLVILVVVG",
#'   S2  = "MTEYKLVVV-LVV"
#' ))
#'
#' # 1) Full span, default reference (last sequence)
#' plotAA(aa)
#'
#' # 2) Choose a reference and zoom
#' plotAA(aa, ref_name = "Ref", xlim = c(3, 8))
#'
#' # 3) Highlight a sequence of interest
#' plotAA(aa, ref_name = "Ref", xlim = c(3, 8),
#'        highlight_seqs = "S1", showLegend = TRUE)
#'
#' @importFrom graphics par plot axis legend rect segments
#' @importFrom grDevices adjustcolor
#' @export
plotAA <- function(aln_aa_set,
                   ref_name = NULL,
                   xlim = NULL,
                   center = NULL, window = NULL,
                   showLegend = FALSE,
                   # 4-group AA color map (same color per class)
                   colorMapAA = c(
                     R = "red2", H = "red2", K = "red2",                         # + charged
                     D = "blue2", E = "blue2",                                  # - charged
                     S = "green3", T = "green3", N = "green3", Q = "green3",      # polar uncharged
                     A = "gold", V = "gold", I = "gold", L = "gold", M = "gold",
                     F = "gold", W = "gold", P = "gold", G = "gold", Y = "gold", C = "gold"  # nonpolar / aromatic / special
                   ),
                   colorFallback = "grey50",
                   highlight_seqs = NULL,
                   highlight_bg = grDevices::adjustcolor("#FFE082", alpha.f = 0.5),
                   lwd_normal = 1,
                   lwd_highlight = 2,
                   seg_half_height = 0.35,
                   y_cex = 0.8,
                   line_col = "grey20",
                   show_only_highlighted = FALSE) {
  
  # Build an AA cache using the same engine
  cache <- compute_cache(
    aln_set       = aln_aa_set,        # AAStringSet allowed with the small tweak above
    ref_name      = ref_name,
    colorMap    = colorMapAA,        # pass AA map into the same arg
    colorFallback = colorFallback
  )
  
  # Draw (full span if no xlim/center+window specified)
  plot_window(
    cache                 = cache,
    xlim                  = xlim,
    center                = center, window = window,
    showLegend            = showLegend,
    highlight_seqs        = highlight_seqs,
    highlight_bg          = highlight_bg,
    lwd_normal            = lwd_normal,
    lwd_highlight         = lwd_highlight,
    seg_half_height       = seg_half_height,
    y_cex                 = y_cex,
    line_col              = line_col,
    show_only_highlighted = show_only_highlighted
  )
  invisible(cache)
}
