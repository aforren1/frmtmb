# Is the distributional fit's vcov() a PERMUTATION of the base matrix,
# or did a number change? Positional comparison cannot tell them apart,
# so match by name and then compare.
TREE <- "C:/Users/adf44/source/r/frmtmb-wt-shapes"
b <- readRDS(file.path(TREE, "dev/shapes-rev-se-base.rds"))
l <- readRDS(file.path(TREE, "dev/shapes-rev-se-lane.rds"))

nmmap <- function(x) x  # brms drops nothing on this fit; names are the map

for (fnm in names(b)) {
  if (fnm == "importance") next
  for (q in c("vcov", "sandwich")) {
    B <- b[[fnm]][[q]]; L <- l[[fnm]][[q]]
    if (inherits(B, "revErr") || inherits(L, "revErr")) next
    cat("\n== ", fnm, " ", q, "\n", sep = "")
    cat("  base rows: ", paste(rownames(B), collapse = " | "), "\n")
    cat("  lane rows: ", paste(rownames(L), collapse = " | "), "\n")
    # try to line the two up by stripping the base's internal decoration
    bn <- rownames(B); ln <- rownames(L)
    cand <- list(ln,
                 sub("^", "", ln),
                 gsub("Intercept", "(Intercept)", ln),
                 gsub("^sigma_", "sigma~", ln))
    hit <- NULL
    for (cc in cand) if (all(cc %in% bn)) { hit <- cc; break }
    if (is.null(hit)) {
      # manual: map brms name -> base name by position-free heuristics
      cat("  [no automatic name map]\n"); next
    }
    Bs <- B[hit, hit, drop = FALSE]
    d <- max(abs(as.numeric(Bs) - as.numeric(L)) /
               pmax(abs(as.numeric(Bs)), 1e-12))
    cat(sprintf("  after name match, max rel diff = %.3e  (%s)\n", d,
                if (identical(as.numeric(Bs), as.numeric(L)))
                  "bitwise identical" else "not bitwise"))
    cat("  base kept rows dropped: ",
        paste(setdiff(bn, hit), collapse = ", "), "\n")
  }
}
