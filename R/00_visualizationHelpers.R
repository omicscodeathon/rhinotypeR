# ===================== Internal topic page (one per file) =====================

#' Internal helpers for SNPeek/plotAA
#'
#' Utilities used by SNPeek/plotAA for caching and fast drawing.
#' These functions are documented for transparency but are **not** part of the
#' public API.
#'
#' @name SNPeek-internal
#' @keywords internal
NULL


# ===================== 1) compareAndColorSequences ============================

#' Compare one reference sequence to others and color substitutions
#'
#' Given an aligned set of sequences (first element treated as the reference),
#' find positions where each non-reference sequence differs and assign colors
#' based on a user-supplied map.
#'
#' @param sequences Character vector (named preferred) of equal-length, aligned
#'   sequences. The first element is treated as the reference.
#' @param colorMap Named character vector mapping symbols (e.g., "A","C","G","T"
#'   or amino-acid one-letter codes) to colors understood by base graphics.
#' @param colorFallback Single color used when a symbol is not present in
#'   \code{colorMap}.
#'
#' @return A list (length = \code{length(sequences) - 1}) of data frames with
#'   columns \code{position} (integer), \code{subsType} (character), and
#'   \code{color} (character), one per non-reference sequence.
#'
#' @seealso \code{\link[=compute_cache]{compute_cache}},
#'   \code{\link[=plot_window]{plot_window}}
#' @examples
#' \dontrun{
#' seqs <- c(Ref="ATGC", S1="ATGT", S2="ATAC")
#' cmap <- c(A="green", T="red", C="blue", G="gold")
#' compareAndColorSequences(seqs, cmap)
#' }
#' 
#' @aliases compareAndColorSequences
#' @rdname SNPeek-internal
#' @keywords internal
compareAndColorSequences <- function(sequences, colorMap, colorFallback = "gray") {
  # Enhanced function to compare sequences and return positions and substitution types
  compareSequences <- function(seqA, seqB) {
    seqAChars <- strsplit(seqA, "")[[1]]
    seqBChars <- strsplit(seqB, "")[[1]]
    differences <- which(seqAChars != seqBChars)
    subsType <- seqBChars[differences]
    data.frame(position = differences, subsType = subsType)
  }
  
  # Compare sequences and assign colors
  diffList <- lapply(2:length(sequences), function(i) {
    diffs <- compareSequences(sequences[[1]], sequences[[i]])
    diffs$color <- colorMap[diffs$subsType]
    diffs$color[is.na(diffs$color)] <- colorFallback
    return(diffs)
  })
  
  return(diffList)
}


# ===================== 2) compute_cache ======================================

#' Build a fast-draw cache for SNPeek/plotAA
#'
#' Precomputes mismatch positions vs. a chosen reference and compact color
#' indices for quick re-drawing of windows without rescanning the alignment.
#'
#' @param aln_set \code{Biostrings::DNAStringSet} or \code{Biostrings::AAStringSet}
#'   containing an aligned, equal-width alignment.
#' @param ref_name Optional character; the sequence name to use as the reference.
#'   Defaults to the last sequence in \code{aln_set}.
#' @param colorMap Named vector mapping symbols (e.g., nucleotides) to colors.
#'   Used to render mismatches. For symbols not present, \code{colorFallback} is used.
#' @param colorFallback Color used when a symbol is not in \code{colorMap}.
#'
#' @return An object of class \code{"SNPeekCache"} (a list) with elements:
#'   \itemize{
#'     \item \code{genome_len} (integer)
#'     \item \code{seq_names} (character)
#'     \item \code{ref_label} (character)
#'     \item \code{ref_row} (integer; row index of the reference)
#'     \item \code{diffs_list} (list of per-sequence mismatch positions & color indices)
#'     \item \code{col_levels} (character; palette including fallback as last element)
#'     \item \code{colorMap}, \code{colorFallback}
#'   }
#'
#' @examples
#' \dontrun{
#' f <- system.file("extdata","test.fasta", package="rhinotypeR")
#' aln <- Biostrings::readDNAStringSet(f)
#' cache <- compute_cache(aln, ref_name = tail(names(aln), 1))
#' }
#'
#' @aliases compute_cache SNPeekCache
#' @rdname SNPeek-internal
#' @keywords internal
compute_cache <- function(aln_set,
                          ref_name = NULL,
                          colorMap,
                          colorFallback = "grey") {
  if (!inherits(aln_set, c("DNAStringSet", "AAStringSet"))) {
    stop("SNPeek/plotAA expects an aligned Biostrings::DNAStringSet or AAStringSet.")
  }
  
  w <- Biostrings::width(aln_set)
  if (length(unique(w)) != 1L) {
    stop("All sequences must be the same width (provide an aligned DNAStringSet).")
  }
  
  seq_names <- names(aln_set)
  genome_len <- w[1]
  
  # Choose reference (default last)
  ref_idx <- if (is.null(ref_name)) length(aln_set) else match(ref_name, seq_names)
  if (is.na(ref_idx)) stop("Reference name not found in alignment: ", ref_name,
                           ". Using the last sequence in the alignment instead.")
  ref_label <- seq_names[ref_idx]
  
  # Reorder so REFERENCE is LAST (stable with existing plotting logic)
  ord_idx <- c(setdiff(seq_along(aln_set), ref_idx), ref_idx)
  aln_ord <- aln_set[ord_idx]
  names_ord <- names(aln_ord)
  ref_row <- length(aln_ord)
  
  # Character matrix (rows = seqs, cols = positions)
  M <- as.matrix(aln_ord)
  
  # Build color index mapping for speed/memory
  nts <- c(names(colorMap))
  pal <- unname(colorMap)
  col_levels <- c(pal, colorFallback)     # last entry is fallback (e.g., N, -)
  
  map_nt_to_col_idx <- function(v) {
    idx <- match(v, nts)
    out <- ifelse(is.na(idx), length(col_levels), idx)  # fallback index = last
    as.integer(out)
  }
  
  # Reference row nucleotides
  ref_nt <- M[ref_row, ]
  
  # Collect mismatch positions and color indices per sequence
  nseq <- nrow(M)
  diffs_list <- vector("list", nseq)
  # placeholder (reference row)
  diffs_list[[ref_row]] <- list(pos = integer(0), col_idx = integer(0))
  
  if (nseq > 1L) {
    for (i in seq_len(nseq - 1L)) {  # all except reference
      neq <- M[i, ] != ref_nt
      if (any(neq)) {
        pos <- which(neq)
        col_idx <- map_nt_to_col_idx(M[i, pos])
        diffs_list[[i]] <- list(pos = pos, col_idx = col_idx)
      } else {
        diffs_list[[i]] <- list(pos = integer(0), col_idx = integer(0))
      }
    }
  }
  
  structure(
    list(
      genome_len = genome_len,
      seq_names  = names_ord,
      ref_label  = ref_label,
      ref_row    = ref_row,
      diffs_list = diffs_list,
      col_levels = col_levels,
      colorMap = colorMap,
      colorFallback = colorFallback
    ),
    class = "SNPeekCache"
  )
}


# ===================== 2b) subset_cache ======================================

#' Subset a precomputed SNPeek cache to selected sequences
#'
#' Keeps only the specified row indices from an \code{SNPeekCache}. Useful for
#' drawing views that focus on a subset (e.g., highlighted sequences).
#'
#' @param cache An object created by \code{\link[=compute_cache]{compute_cache}}
#'   (class \code{"SNPeekCache"}).
#' @param keep_idx Integer vector of 1-based row indices to retain. Values
#'   outside the valid range are ignored.
#'
#' @return A pruned \code{"SNPeekCache"} object containing only the selected
#'   sequences. The original metadata (e.g., \code{ref_row}) is preserved and
#'   is not used directly by the drawing code.
#'
#' @examples
#' \dontrun{
#' cache2 <- subset_cache(cache, c(1, 3, 5))
#' }
#' 
#' @aliases subset_cache
#' @rdname SNPeek-internal
#' @keywords internal
subset_cache <- function(cache, keep_idx) {
  stopifnot(inherits(cache, "SNPeekCache"))
  keep_idx <- sort(unique(keep_idx))
  keep_idx <- keep_idx[keep_idx >= 1 & keep_idx <= length(cache$seq_names)]
  cache2 <- cache
  cache2$seq_names  <- cache$seq_names[keep_idx]
  cache2$diffs_list <- cache$diffs_list[keep_idx]
  # ref_row/meta can stay pointing to original ref; not used directly in drawing
  cache2
}


# ===================== 3) plot_window ========================================

#' Draw a windowed mismatch view from an SNPeek cache
#'
#' Low-level renderer used by higher-level plotting helpers. Plots mismatch
#' ticks relative to a cached reference over a specified genomic span, with
#' optional highlighting and a legend.
#'
#' @param cache An object from \code{\link[=compute_cache]{compute_cache}}
#'   (class \code{"SNPeekCache"}).
#' @param xlim Integer length-2 vector giving the genomic range to display.
#'   If omitted, you can provide \code{center} and \code{window} for quick zoom.
#' @param center Integer center position for quick zoom (used with \code{window}).
#' @param window Integer window width for quick zoom (used with \code{center}).
#' @param showLegend Logical; draw a small legend for the color map.
#' @param highlight_seqs Character vector of sequence names to visually
#'   highlight (background shading and thicker ticks).
#' @param highlight_bg Background color for highlighted rows (can include alpha).
#' @param lwd_normal Numeric line width for non-highlighted ticks.
#' @param lwd_highlight Numeric line width for highlighted ticks.
#' @param seg_half_height Half-height of mismatch ticks (in plot y-units).
#' @param y_cex Character expansion for y-axis (sequence name) labels.
#' @param line_col Color for axes/frame elements.
#' @param show_only_highlighted Logical; if \code{TRUE} and \code{highlight_seqs}
#'   are provided, only those rows are drawn.
#'
#' @return Invisibly returns a list with \code{xlim}, \code{n_sequences}, and
#'   \code{highlighted} sequence names.
#'
#' @seealso \code{\link[=compute_cache]{compute_cache}}
#' @examples
#' \dontrun{
#' cache <- compute_cache(aln)
#' plot_window(cache, center = 200, window = 150,
#'             highlight_seqs = cache$seq_names[1:3], showLegend = TRUE)
#' }
#'
#' @aliases plot_window
#' @rdname SNPeek-internal
#' @keywords internal
plot_window <- function(cache,
                        xlim = NULL,
                        center = NULL, window = NULL,
                        showLegend = FALSE,
                        highlight_seqs = NULL,
                        highlight_bg = "#FFE08280",   # translucent amber
                        lwd_normal = 1,
                        lwd_highlight = 2,
                        seg_half_height = 0.35,
                        y_cex = 0.8,
                        line_col = "grey20",
                        show_only_highlighted = FALSE) {
  stopifnot(inherits(cache, "SNPeekCache"))
  
  genome_len <- cache$genome_len
  if (!is.null(center) && !is.null(window)) {
    half <- floor(window / 2)
    xlim <- c(max(1L, center - half), min(genome_len, center + half))
  }
  if (is.null(xlim)) xlim <- c(1L, genome_len)
  
  # Map highlight sequence names to row indices
  hi_idx <- integer(0)
  if (!is.null(highlight_seqs)) {
    hi_idx <- match(highlight_seqs, cache$seq_names)
    if (any(is.na(hi_idx))) {
      warning("Some 'highlight_seqs' not found and will be ignored: ",
              paste(highlight_seqs[is.na(hi_idx)], collapse = ", "))
      hi_idx <- hi_idx[!is.na(hi_idx)]
    }
  }
  
  # Optionally reduce to only highlighted rows (still mismatches vs cached ref)
  if (show_only_highlighted && length(hi_idx)) {
    cache <- subset_cache(cache, hi_idx)
    # Renumber hi_idx in this reduced view:
    hi_idx <- seq_along(cache$seq_names)
  } else if (show_only_highlighted && !length(hi_idx)) {
    warning("'show_only_highlighted = TRUE' but no valid 'highlight_seqs' given; showing all.")
  }
  
  oldPar <- par(mar = c(5, 8, 4, 2) + 0.1)
  on.exit(par(oldPar), add = TRUE)
  
  plot(NULL,
       xlim = xlim,
       ylim = c(0.5, length(cache$seq_names) + 0.5),
       type = "n",
       xlab = paste0("Genome position of ", cache$ref_label, " (reference)"),
       ylab = "", yaxt = "n")
  
  axis(2, at = seq_along(cache$seq_names), labels = cache$seq_names, las = 2, cex.axis = y_cex)
  
  # Soft background bands for highlighted rows (draw first)
  if (length(hi_idx)) {
    usr <- par("usr")  # (x1, x2, y1, y2)
    for (i in hi_idx) {
      rect(xleft = usr[1], ybottom = i - 0.5,
           xright = usr[2], ytop = i + 0.5,
           col = highlight_bg, border = NA)
    }
  }
  
  # Draw diffs, subsetting to xlim first
  xmin <- xlim[1]; xmax <- xlim[2]
  for (i in seq_along(cache$diffs_list)) {
    d <- cache$diffs_list[[i]]
    if (length(d$pos)) {
      keep <- (d$pos >= xmin) & (d$pos <= xmax)
      if (any(keep)) {
        px   <- d$pos[keep]
        cols <- cache$col_levels[d$col_idx[keep]]
        lw   <- if (i %in% hi_idx) lwd_highlight else lwd_normal
        segments(x0 = px, y0 = i - seg_half_height,
                 x1 = px, y1 = i + seg_half_height,
                 col = cols, lwd = lw)
      }
    }
  }
  
  if (showLegend) {
    leg_labs <- c(names(cache$colorMap), "Other")
    leg_cols <- c(unname(cache$colorMap), cache$colorFallback)
    legend("topleft", inset = c(0.8, 0),
           legend = leg_labs,
           fill   = leg_cols,
           cex = 0.6, bty = "n", box.col = "gray",
           bg = "white")
  }
  
  invisible(list(
    xlim = xlim,
    n_sequences = length(cache$seq_names),
    highlighted = cache$seq_names[hi_idx]
  ))
}


