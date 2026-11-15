# Lane predraw2: how many test draws the check needs. For each fit, the
# excess e of draw_logratio() at 1000 draws of the fit's own law (seed 7,
# not the check's), then, for n0 = 10, 20 and 40, the share of 2000 random
# n0-subsets whose median is above draw_check_tol: the chance that a check
# with n0 test draws fires on that fit.
# Usage: Rscript dev/predraw2-ncheck.R <lib> <group> <out.tsv>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
source("dev/predraw2-cases.R")
group <- args[2]
out <- args[3]
ids <- if (group == "mixg40") {
  paste0("mixgb:", 1:40)
} else if (group == "S") {
  c(paste0("S0.1:", 1:25), paste0("S0.2:", 1:25))
} else {
  case_ids(group)
}
fit_of <- function(id) {
  if (startsWith(id, "mixgb:")) {
    d <- gen("mixture_gaussian", "boundary", as.integer(sub(".*:", "", id)))
    return(fq(frm(bf(y ~ x), family = mixture(gaussian(), gaussian()),
                  data = d)))
  }
  case_fit(id)$fit
}
rows <- list()
for (id in ids) {
  fit <- tryCatch(fit_of(id), error = function(e) NULL)
  if (is.null(fit)) next
  # the law predict() builds: run once, then rebuild the pieces
  set.seed(7)
  dr <- fq(ns$predict_par_drawer(fit, TRUE, 2L))
  ds <- ns$fit_draw_space(fit)
  V <- ds$V
  if (is.null(V) || !all(is.finite(V))) next
  L <- tryCatch(chol(V), error = function(e) {
    ev <- eigen(V, symmetric = TRUE)
    sqrt(pmax(ev$values, 0)) * t(ev$vectors)
  })
  if (is.null(L)) next
  law <- fit$cache[["draw_edge_law"]]
  if (!is.list(law)) law <- NULL
  v0 <- ns$fit_outer_vector(fit, ds$map)
  P <- length(v0)
  dw <- list(L = L, law = law, v0 = v0, held = integer(0),
             off = P - length(fit$opt$par))
  set.seed(7)
  Z <- matrix(stats::rnorm(P * 1000), P, 1000)
  D <- v0 + if (is.null(law)) crossprod(L, Z) else
    ns$draw_edge_apply(law, L, Z)
  r <- ns$draw_logratio(fit, dw, D, Z)
  e <- r$e
  set.seed(11)
  pf <- vapply(c(10L, 20L, 40L), function(n0) {
    mean(replicate(2000, stats::median(sample(e, n0)) > ns$draw_check_tol))
  }, 0)
  # sequential: 20 draws, and 20 more where their median is within a
  # factor of 2 of the threshold, deciding on all 40
  tol <- ns$draw_check_tol
  ps <- mean(replicate(2000, {
    s <- sample(e, 40L)
    m20 <- stats::median(s[1:20])
    if (m20 >= tol / 2 && m20 <= 2 * tol) stats::median(s) > tol else
      m20 > tol
  }))
  rows[[id]] <- data.frame(id = id, med1000 = stats::median(e),
                           p10 = pf[1], p20 = pf[2], p40 = pf[3], pseq = ps)
  cat(sprintf(paste("%-14s med1000 %8.3f  P(fire) n10 %.4f n20 %.4f n40",
                    "%.4f seq20+20 %.4f\n"),
              id, stats::median(e), pf[1], pf[2], pf[3], ps))
}
res <- do.call(rbind, rows)
utils::write.table(res, out, sep = "\t", row.names = FALSE, quote = FALSE)
