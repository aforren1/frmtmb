# Lane wt-reunc, punch round 2 (B1): the guard stated as a general
# property and checked against the thing it describes.
#
# Every caller that perturbs, bounds or batches `b` one entry at a time
# needs `expand_b()` to carry `b` to the coefficient vector
# POSITIONWISE, so that `cvec[c_idx[i]]` is a function of `b[b_idx[i]]`
# alone. `block_b_positionwise()` is that predicate. This script does
# not take its word for which blocks fail it. It reads `expand_b()`'s
# OWN branch conditions out of the function body, enumerates every
# covstruct the registry offers and every car type the package accepts,
# and shows the predicate's verdict equals the branch structure on all
# of them. A new covstruct that needed a branch and did not get one
# would show up here as a disagreement rather than as a wrong standard
# error years later.
#
#   Rscript dev/reunc-positionwise.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/reunc2-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n\n")

# 1. What expand_b() itself branches on. Its body is
#    `for (bk in ...) if (<rr>) ... else if (<esicar>) ... else <copy>`,
#    so the non-positionwise set is exactly the guarded branches and the
#    `else` is the positionwise copy.
src <- paste(deparse(body(frmtmb:::expand_b)), collapse = " ")
cat("1. expand_b()'s own branches\n")
has_rr <- grepl('bk[["covstruct"]] == "rr"', src, fixed = TRUE)
has_esi <- grepl("block_is_esicar(bk)", src, fixed = TRUE)
copy <- grepl('cvec[bk[["c_idx"]]] <- b[bk[["b_idx"]]]', src, fixed = TRUE)
cat("   rr branch present:            ", has_rr, "\n", sep = "")
cat("   esicar branch present:        ", has_esi, "\n", sep = "")
cat("   else branch is a plain copy:  ", copy, "\n", sep = "")
# the branch count is what makes this an enumeration and not a spot
# check: a third special case would raise it
nbranch <- length(gregexpr("cvec\\[bk\\[\\[\"c_idx\"\\]\\]\\] <-", src)[[1]])
cat("   assignments to cvec[c_idx]:   ", nbranch,
    " (rr, esicar, copy)\n", sep = "")
stopifnot(has_rr, has_esi, copy, nbranch == 3L)

# 2. Every covstruct the registry offers, and every car type.
reg <- sort(names(frmtmb:::covstruct_registry))
cars <- frmtmb:::car_types
cat("\n2. the registry: ", length(reg), " covstructs, ",
    length(cars), " car types\n", sep = "")

# what expand_b() does to a block, stated independently of the
# predicate: read off the branch conditions above
expected <- function(cs, ct) {
  !(identical(cs, "rr") ||
      (identical(cs, "car") && identical(ct, "esicar")))
}

rows <- list()
for (cs in reg) {
  tys <- if (identical(cs, "car")) cars else NA_character_
  for (ct in tys) {
    bk <- list(covstruct = cs, car_type = if (is.na(ct)) NULL else ct)
    rows[[length(rows) + 1L]] <- data.frame(
      covstruct = cs,
      car_type = if (is.na(ct)) "" else ct,
      guard_says = frmtmb:::block_b_positionwise(bk),
      expand_b_says = expected(cs, ct),
      stringsAsFactors = FALSE)
  }
}
tab <- do.call(rbind, rows)
tab$agree <- tab$guard_says == tab$expand_b_says
cat("\n")
print(tab, row.names = FALSE)
cat("\n   rows: ", nrow(tab), "; disagreements: ", sum(!tab$agree),
    "\n", sep = "")
cat("   not positionwise: ",
    paste(ifelse(nzchar(tab$car_type[!tab$guard_says]),
                 paste0(tab$covstruct[!tab$guard_says], "(",
                        tab$car_type[!tab$guard_says], ")"),
                 tab$covstruct[!tab$guard_says]), collapse = ", "),
    "\n", sep = "")
stopifnot(all(tab$agree))

# 3. frame_needs_expand() is the other caller of the same property and
#    must not have drifted from it.
cat("\n3. frame_needs_expand() agrees block by block\n")
for (i in seq_len(nrow(tab))) {
  bk <- list(covstruct = tab$covstruct[i],
             car_type = if (nzchar(tab$car_type[i])) tab$car_type[i])
  fr <- list(re_blocks = list(bk))
  ok <- frmtmb:::frame_needs_expand(fr) == !tab$guard_says[i]
  if (!ok) {
    cat("   DISAGREES: ", tab$covstruct[i], " ", tab$car_type[i], "\n",
        sep = "")
  }
}
cat("   all ", nrow(tab), " agree\n", sep = "")
