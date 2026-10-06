# Build the "before" arms of the 0.68.0 merge's resolutions: copies of
# core with one resolution undone, each installed into a library of its
# own, so that every test written for a resolution is seen to fail
# before it and pass after it (dev/round-20261005.md, "Conflicts").
#
#   Rscript dev/rel068-mutants.R <out-dir> <conflicted priors.R>
#
# Arms (each from the release tree as it stands, then one change):
#   naive   the textual merge's own sides: R/priors.R hunks 1 and 2 from
#           lane fixes, hunk 3 from lane ordmix (no coef_k); the compat
#           cells, ord_thres_linpred() and ord_internal_labels() as the
#           lanes left them
#   h3fix   R/priors.R hunk 3 from lane fixes (offset without dpar)
#   h3mix   R/priors.R hunk 3 from lane ordmix (no coef_k)
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
conflicted <- args[2]
root <- "C:/Users/adf44/source/r/frmtmb-wt-release"
base <- "9e902909"
side <- function(lines, picks) {
  # picks: one of "ours" / "theirs" per conflict hunk, in order
  res <- character(0)
  h <- 0L
  mode <- "out"
  for (l in lines) {
    if (startsWith(l, "<<<<<<< ")) {
      h <- h + 1L
      mode <- "ours"
    } else if (startsWith(l, "=======") && mode == "ours") {
      mode <- "theirs"
    } else if (startsWith(l, ">>>>>>> ")) {
      mode <- "out"
    } else if (mode == "out" || mode == picks[h]) {
      res <- c(res, l)
    }
  }
  res
}
cur <- function(f) readLines(file.path(root, f), warn = FALSE)
from_lane <- function(lane, f) {
  readLines(file.path("C:/Users/adf44/source/r", paste0("frmtmb-wt-", lane),
                      f), warn = FALSE)
}
# one function's definition swapped for another file's
swap_fun <- function(lines, donor, name) {
  rng <- function(x) {
    s <- grep(paste0("^", name, " <- function"), x)
    stopifnot(length(s) == 1L)
    e <- s
    while (!identical(x[e], "}")) e <- e + 1L
    c(s, e)
  }
  a <- rng(lines)
  b <- rng(donor)
  c(lines[seq_len(a[1] - 1L)], donor[b[1]:b[2]],
    lines[(a[2] + 1L):length(lines)])
}
hunk3 <- function(lines, which) {
  s <- grep("offset = if (!ob$shared) {", lines, fixed = TRUE)
  stopifnot(length(s) == 1L)
  e <- s + 4L
  stopifnot(grepl("coef_k = k)))", lines[e], fixed = TRUE))
  new <- switch(which,
    fixes = c("                       offset = ordinal_center_offset(frame, rspec),",
              "                       lb = s$lb, ub = s$ub, coef_k = k)))"),
    ordmix = c(lines[s:(e - 1L)],
               "                       lb = s$lb, ub = s$ub)))"))
  c(lines[seq_len(s - 1L)], new, lines[(e + 1L):length(lines)])
}
mk <- function(arm, edit) {
  d <- file.path(out, arm)
  unlink(d, recursive = TRUE)
  dir.create(d, recursive = TRUE)
  for (x in c("DESCRIPTION", "NAMESPACE", "R", "inst", "data")) {
    file.copy(file.path(root, x), d, recursive = TRUE)
  }
  edit(d)
  cat("built", arm, "\n")
}
wr <- function(d, f, lines) writeLines(lines, file.path(d, f))
mk("naive", function(d) {
  wr(d, "R/priors.R", side(readLines(conflicted, warn = FALSE),
                            c("ours", "ours", "theirs")))
  # compat.R: the merge took both lanes' rows; only the three cells
  # differ, so take the release file and put the lane's loop back
  rel <- cur("R/compat.R")
  rm <- grep('^  r\\("(disc|equidistant|sum_to_zero)", "mixture"', rel)
  rel <- rel[-c(rm, rm + 1L)]
  at <- grep('for (m in c("fitted", "predict", "simulate", "residuals_osa"))',
             rel, fixed = TRUE)[1]
  rel <- append(rel, c(
    '    r(st, "mixture", "refused",',
    '      "Refused with every ordinal family: an ordinal family is not a mixture component.")'),
    after = at - 1L)
  wr(d, "R/compat.R", rel)
  wr(d, "R/predict.R", swap_fun(cur("R/predict.R"),
                                from_lane("fixes", "R/predict.R"),
                                "ord_thres_linpred"))
  wr(d, "R/confint.R", swap_fun(cur("R/confint.R"),
                                from_lane("fixes", "R/confint.R"),
                                "ord_internal_labels"))
})
mk("h3fix", function(d) wr(d, "R/priors.R", hunk3(cur("R/priors.R"),
                                                  "fixes")))
mk("h3mix", function(d) wr(d, "R/priors.R", hunk3(cur("R/priors.R"),
                                                  "ordmix")))
cat("MUTANTS DONE\n")
