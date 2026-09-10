# rev-gddm: how much of the tolerance the accepted designs actually
# need, measured on the lane's own sweep.
#
# frmtmb_register_frame_check() is public, so a recorder can be added
# AFTER gd_check_condition_constancy(). Checks run in registration
# order and the gddm one throws, so the recorder sees only frames the
# guard ACCEPTED. For every numeric model-frame column of a gddm
# response it records
#
#     max |row - its condition's first row| / max |finite entries|
#
# which is the quantity the shipped constant 1e-8 is compared against.
# A superset of the columns the check itself looks at, so the number
# here bounds the one that matters.
#
# Then it reruns dev/gddm-scripts/gddm-sweep.R unchanged (seed 202609)
# so the false-alarm and drop counts are reproduced in the same
# process.

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/rev-gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

REC <- new.env(parent = emptyenv())
REC$rel <- numeric(0)
REC$where <- character(0)
REC$frames <- 0L

recorder <- function(spec, frame) {
  sp <- frame[["spec"]] %||% spec
  mf <- frame[["data_frame"]]
  av <- frame[["aterm_values"]]
  if (is.null(sp) || is.null(mf) || !nrow(mf)) return(invisible(NULL))
  gd <- vapply(sp[["responses"]],
               function(r) !is.null(r[["family"]][["gddm"]]), TRUE)
  if (!any(gd)) return(invisible(NULL))
  REC$frames <- REC$frames + 1L
  for (rn in names(sp[["responses"]])[gd]) {
    av_r <- av[[rn]]
    cnd <- if (!is.null(av_r[["dec"]])) av_r[["vint1"]] else
      av_r[["vint2"]]
    if (is.null(cnd) || length(cnd) != nrow(mf)) next
    lev <- sort(unique(cnd)); gi <- match(cnd, lev)
    first <- match(seq_along(lev), gi)
    # the same restriction the check makes: only columns built from a
    # variable some dpar's formula names
    resp <- sp[["responses"]][[rn]]
    vars <- unique(unlist(lapply(names(resp[["dpars"]]), function(dn) {
      dp <- resp[["dpars"]][[dn]]
      v <- all.vars(dp[["rhs"]])
      nl <- dp[["nl_body"]]
      if (!is.null(nl)) v <- c(v, setdiff(all.vars(nl),
                                          resp[["nlpars"]]))
      v
    })))
    if (!length(vars)) next
    cv <- lapply(names(mf), function(cn)
      tryCatch(all.vars(str2lang(cn)), error = function(e) cn))
    for (k in seq_along(mf)) {
      if (!length(intersect(cv[[k]], vars))) next
      v <- mf[[k]]
      if (is.data.frame(v) || is.matrix(v)) v <- as.matrix(v) else
        if (!is.numeric(v)) next
      m <- as.matrix(v)
      if (!is.numeric(m)) next
      r <- m[first[gi], , drop = FALSE]
      for (j in seq_len(ncol(m))) {
        fin <- is.finite(m[, j]) & is.finite(r[, j])
        if (!any(fin)) next
        sc <- max(abs(m[fin, j]))
        if (sc <= 0) next
        rel <- max(abs(m[fin, j] - r[fin, j])) / sc
        REC$rel <- c(REC$rel, rel)
        REC$where <- c(REC$where, paste0(names(mf)[k], "[", j, "]"))
      }
    }
  }
  invisible(NULL)
}
frmtmb_register_frame_check(recorder)

`%||%` <- function(x, y) if (is.null(x)) y else x
source("C:/Users/adf44/source/r/frmtmb-wt-gddm/dev/gddm-scripts/gddm-sweep.R",
       local = new.env())

cat("\n\n== recorder: columns of ACCEPTED gddm frames\n")
cat("  gddm frames accepted:", REC$frames, "\n")
cat("  column-condition comparisons recorded:", length(REC$rel), "\n")
nz <- REC$rel[REC$rel > 0]
cat("  of which exactly zero:", sum(REC$rel == 0),
    "  nonzero:", length(nz), "\n")
if (length(nz)) {
  o <- order(nz, decreasing = TRUE)
  cat("  the nonzero relative deviations, largest first:\n")
  for (i in head(o, 12L)) {
    cat(sprintf("    %-24s %.3e\n", REC$where[REC$rel > 0][i], nz[i]))
  }
  cat(sprintf("  worst: %.3e   shipped tolerance: 1.0e-08",
              max(nz)))
  cat(sprintf("   headroom %.3g\n", 1e-8 / max(nz)))
  for (cst in c(1e-9, 1e-10, 1e-11, 1e-12, 1e-13)) {
    cat(sprintf("    a constant of %.0e would false-alarm on %d of %d",
                cst, sum(nz > cst), length(REC$rel)))
    cat(" accepted comparisons\n")
  }
}
