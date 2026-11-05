# The proposed fix as a mutant, to show the step is the cause: each beta
# coordinate's step is 1e-4 of a unit of the predictor, from the largest absolute
# value of its design column, and not 1e-4 of the coefficient's own units. Copy of
# the source under rev-impsamp-src/fix_x, nothing in the lane's tree changes.
WT <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp"
dst <- "C:/Users/adf44/source/r/rev-impsamp-src/fix_x"
unlink(dst, recursive = TRUE)
dir.create(file.path(dst, "R"), recursive = TRUE)
dir.create("C:/Users/adf44/source/r/rev-impsamp-lib/fix_x", recursive = TRUE,
           showWarnings = FALSE)
for (f in c("DESCRIPTION", "NAMESPACE")) file.copy(file.path(WT, f), dst)
file.copy(list.files(file.path(WT, "R"), full.names = TRUE), file.path(dst, "R"))
for (d in c("inst", "data", "man")) file.copy(file.path(WT, d), dst, recursive = TRUE)
p <- file.path(dst, "R", "importance.R")
s <- readLines(p)
i <- grep("sl <- imp_anchor_slopes(lap_obj, full, lay, hlik, prec)", s, fixed = TRUE)
stopifnot(length(i) == 1L)
s <- append(s, c(
  "  attr(lay, \"xscale\") <- unlist(lapply(frame[[\"linpreds\"]], function(lp) {",
  "    if (is.null(lp[[\"X\"]]) || !ncol(lp[[\"X\"]])) return(NULL)",
  "    apply(abs(as.matrix(lp[[\"X\"]])), 2L, max)",
  "  }))"), i - 1L)
i <- grep("for (j in seq_along(fix)) {", s, fixed = TRUE)
stopifnot(length(i) == 1L)
s <- append(s, c(
  "  xs <- attr(lay, \"xscale\")",
  "  ib <- which(names(full)[fix] == \"beta\")",
  "  hbeta <- rep(NA_real_, length(fix))",
  "  if (length(xs) == length(ib)) hbeta[ib] <- imp_cross_step / pmax(xs, 1e-8)"),
  i - 1L)
i <- grep("h <- imp_cross_step * max(1, abs(full[[fix[[j]]]]))", s, fixed = TRUE)
stopifnot(length(i) == 1L)
s <- append(s, "    if (is.finite(hbeta[[j]])) h <- hbeta[[j]]", i)
writeLines(s, p)
cat("edited\n")
