# Lane optstart, item R1: the fits whose standard errors come from the
# finite-difference Hessian in natural units (par_units other than 1),
# from the HESS lines dev/optstart-profile.R writes. Per run: how many,
# how many whose route Hessian is not positive definite where a
# Richardson Jacobian of the gradient is, and the distribution of the
# largest relative SE error against that Jacobian. Matched by test file,
# key and occurrence between two runs; the worst fits are listed.
#   Rscript dev/optstart-hess-sum.R <base run> <lane run>
a <- commandArgs(TRUE)
rd <- function(run) {
  fs <- Sys.glob(file.path("dev/optstart-log", paste0("prof-", run), "*.txt"))
  # each HESS line follows the FIT line of its fit in the same file;
  # carry that fit's convergence code and warnings onto it
  L <- unlist(lapply(fs, function(f) {
    z <- readLines(f, warn = FALSE)
    last <- ""
    out <- character(0)
    for (l in z) {
      if (startsWith(l, "FIT\t")) last <- l
      if (startsWith(l, "HESS\t")) {
        q <- strsplit(last, "\t", fixed = TRUE)[[1]]
        out <- c(out, paste(l, if (length(q) >= 9) q[6] else "NA",
                            if (length(q) >= 9) q[9] else "-", sep = "\t"))
      }
    }
    out
  }))
  L <- L[startsWith(L, "HESS\t")]
  p <- strsplit(L, "\t", fixed = TRUE)
  X <- data.frame(tag = vapply(p, `[`, "", 2), key = vapply(p, `[`, "", 3),
                  npar = as.integer(vapply(p, `[`, "", 4)),
                  tau = as.numeric(vapply(p, `[`, "", 5)),
                  err = suppressWarnings(as.numeric(vapply(p, `[`, "", 6))),
                  pd_route = vapply(p, `[`, "", 7) == "TRUE",
                  pd_rich = vapply(p, `[`, "", 8) == "TRUE",
                  nbad = as.integer(vapply(p, `[`, "", 9)),
                  code = vapply(p, `[`, "", 10),
                  warned = vapply(p, `[`, "", 11) != "-",
                  stringsAsFactors = FALSE)
  X$occ <- stats::ave(seq_len(nrow(X)), X$tag, X$key, FUN = seq_along)
  X
}
summ <- function(X, lab) {
  ok <- X$pd_rich
  q <- ok & !X$warned & X$code == "0"
  cat(sprintf(paste0("%s: %d fits (%d test files) | Richardson positive ",
                     "definite on %d, %d of them at code 0 without a ",
                     "warning | of those, route not positive definite %d, ",
                     "SE error > 1e-3 %d, > 1e-2 %d, > 0.1 %d | largest ",
                     "rise tau > 1e-2 on %d\n"),
              lab, nrow(X), length(unique(X$tag)), sum(ok), sum(q),
              sum(q & !X$pd_route),
              sum(q & X$pd_route & !is.na(X$err) & X$err > 1e-3),
              sum(q & X$pd_route & !is.na(X$err) & X$err > 1e-2),
              sum(q & X$pd_route & !is.na(X$err) & X$err > 0.1),
              sum(X$tau > 1e-2, na.rm = TRUE)))
}
B <- rd(a[1])
N <- rd(a[2])
summ(B, a[1])
summ(N, a[2])
M <- merge(B, N, by = c("tag", "key", "occ"), suffixes = c(".b", ".n"))
cat("matched:", nrow(M), "| base only:", nrow(B) - nrow(M), "| lane only:",
    nrow(N) - nrow(M), "\n")
off <- function(pd, err) !pd | (!is.na(err) & err > 1e-3)
bad <- M[M$pd_rich.b & (off(M$pd_route.b, M$err.b) |
                        off(M$pd_route.n, M$err.n)), ]
bad <- bad[order(-pmax(bad$err.b, bad$err.n, na.rm = TRUE)), ]
cat("\nfits where either build's route is not positive definite or off by",
    "more than 1e-3 (Richardson positive definite):", nrow(bad), "\n")
for (i in seq_len(nrow(bad))) {
  cat(sprintf(paste0("  %-40s n %2d | base pd %s err %.3g tau %.3g code %s",
                     "%s | lane pd %s err %.3g tau %.3g code %s%s\n"),
              bad$tag[i], bad$npar.b[i], bad$pd_route.b[i], bad$err.b[i],
              bad$tau.b[i], bad$code.b[i], if (bad$warned.b[i]) " w" else "",
              bad$pd_route.n[i], bad$err.n[i], bad$tau.n[i], bad$code.n[i],
              if (bad$warned.n[i]) " w" else ""))
}
