# Author: Kyun-Seop Bae k@acr.kr
# Internal plotting helpers. Their names start with a dot, which the
# exportPattern in NAMESPACE does not match, so they stay out of the public API.

.pkPalette = function(n)
{
# Distinguishable colours for an arbitrary number of series.
# The first 18 are the palette pkr has always used, so existing studies look
# unchanged; beyond that a generated qualitative palette takes over. Indexing
# past the end of a fixed vector yields NA, and an NA colour draws nothing at
# all, which used to make every series after the 18th silently disappear.
  base18 = c("navy", "orangered2", "yellowgreen", "magenta", "turquoise4", "tomato4",
             "mediumblue", "olivedrab4", "yellow3", "indianred", "lightcoral",
             "seashell3", "hotpink3", "midnightblue", "peru", "plum2", "khaki3", "lightgreen")
  if (!is.finite(n) || n < 1) return(character(0))
  if (n <= length(base18)) return(base18[seq_len(n)])
  return(grDevices::hcl.colors(n, palette="Dark 3"))
}

.pkPch = function(n)
{
# Symbols cycle more slowly than colours, so two series sharing a hue in a
# large study can still be told apart.
  if (!is.finite(n) || n < 1) return(integer(0))
  return(rep(c(16, 17, 15, 18, 1, 2, 0, 5), length.out=n))
}

.pkLogAxis = function(lo, hi, maxTicks=9)
{
# Decade-aligned semi-log axis from a concentration range on the original scale.
# Returns ylim (log10 units) enclosing the ticks, so at least one labelled tick
# is always inside the plot. The previous trunc()-based code could place its
# only tick outside ylim and leave the y axis blank.
  lo = suppressWarnings(log10(lo))
  hi = suppressWarnings(log10(hi))
  if (!is.finite(lo) || !is.finite(hi)) return(NULL)
  if (hi < lo) { tmp = lo; lo = hi; hi = tmp }
  l = floor(lo)
  u = ceiling(hi)
  if (u <= l) u = l + 1
  by = max(1, ceiling((u - l)/maxTicks))
  return(list(ylim=c(l, u), at=seq(l, u, by=by)))
}

.pkDrawLogAxis = function(ax)
{
  if (is.null(ax)) return(invisible(NULL))
  axis(2, at=ax$at, labels=as.expression(lapply(ax$at, function(i) bquote(10^.(i)))))
}

.pkSafeName = function(x, default="Data")
{
# File names are built from deparse(substitute(concData)), so an expression such
# as d[d$ID != "S001", ] would otherwise produce a path the OS rejects.
  x = paste(as.character(x), collapse=" ")
  x = gsub("[[:cntrl:]]", " ", x)
  x = gsub("[\\\\/:*?\"<>|]", "_", x)
  x = gsub("[[:space:]]+", " ", x)
  x = Trim(x)
  x = sub("[. ]+$", "", x)             # Windows silently drops these at the end
  if (nchar(x) > 60) x = Trim(substr(x, 1, 60))
  if (!nzchar(x)) x = default
  return(x)
}

.pkGrid = function(n, maxCol=4)
{
# Panel grid for n treatments. A single row of n panels made both the margin
# request (mar grew with n) and the raster width unbounded.
  if (!is.finite(n) || n < 1) return(c(1, 1))
  nc = min(n, maxCol)
  return(c(ceiling(n/nc), nc))
}

.pkRaster = function(file, nr, nc, panel=5, legend=0, res=150, maxIn=24)
{
# Open a raster device sized from the panel grid, in inches rather than pixels,
# and capped so that many treatments cannot ask for a raster of tens of
# thousands of pixels a side. Falls back to PNG where the build has no TIFF.
  w = nc*panel + legend
  h = nr*panel
  s = min(1, maxIn/w, maxIn/h)
  w = w*s
  h = h*s
  if (isTRUE(capabilities("tiff"))) {
    tiff(filename=file, width=w, height=h, units="in", res=res, compression="lzw")
  } else {
    png(filename=sub("[.]tiff?$", ".png", file), width=w, height=h, units="in", res=res)
  }
  return(invisible(NULL))
}

.pkLegendPanel = function(labels, col, pch=16, lty=1, title="Subject ID",
                          maxCex=0.85, minCex=0.35)
{
# Draw a series legend in its own layout cell, adding columns and shrinking the
# text until it fits. Because the cell is not the plot region, the legend can no
# longer run over the curves or off the figure. When even the smallest setting
# will not fit, report the count rather than draw an unreadable legend.
  op = par(mar=c(0, 0, 0, 0))
  on.exit(par(op))
  plot.new()
  n = length(labels)
  if (n == 0) return(invisible(FALSE))
  widest = labels[which.max(nchar(labels))]
  for (cx in seq(maxCex, minCex, by=-0.05)) {
    for (k in seq_len(4)) {
      rows = ceiling(n/k)
      hNeed = (rows + 2) * strheight("Mg", cex=cx) * 1.35
      wNeed = k * (strwidth(widest, cex=cx) + strwidth("MMM", cex=cx))
      if (hNeed <= 1 & wNeed <= 1) {
        legend("center", legend=labels, col=col, pch=pch, lty=lty, ncol=k,
               cex=cx, title=title, bty="n", xpd=NA)
        return(invisible(TRUE))
      }
    }
  }
  text(0.5, 0.5, paste0(title, "\n\n", n, " series\n(too many to label)"), cex=maxCex)
  return(invisible(FALSE))
}

.pkNominalTime = function(tv, nbins)
{
# Nominal sampling grid for the mean profile. bins.greedy() can fail - more bins
# than distinct times, for one - so fall back to the distinct sampling times,
# which is what a nominal schedule usually is anyway.
  ut = sort(unique(tv[is.finite(tv)]))
  nbins = min(nbins, length(ut))
  b = NULL
  if (nbins >= 2 && requireNamespace("binr", quietly=TRUE)) {
    b = tryCatch(binr::bins.greedy(tv, nbins=nbins, naive=TRUE), error=function(e) NULL)
  }
  if (is.null(b) || is.null(b$binlo) || length(b$binlo) == 0) {
    ct = as.vector(table(tv)[as.character(ut)])
    return(list(t=ut, ct=ct))
  }
  lo = as.vector(b$binlo)
  return(list(t=ifelse(lo < 5, round(lo, 1), round(lo)), ct=as.vector(b$binct)))
}
