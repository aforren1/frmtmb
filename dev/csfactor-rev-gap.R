# Why the rank test misses mo(m) + cs(m): the stored X column for mo() is
# a zero placeholder, so rnk(cbind(1, X)) does not count it and the cs
# column that mo() can reproduce still raises the rank.
# Also: the compat registry row, and the importance refusal.
#   Rscript dev/csfactor-rev-gap.R <lib> <tag>
args <- commandArgs(TRUE)
LIB <- args[[1L]]; TAG <- args[[2L]]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== build", TAG, ":", as.character(packageVersion("frmtmb")), "\n")
short <- function(e) substr(gsub("[\r\n]+", " ", conditionMessage(e)),
                           1L, 400L)

set.seed(1907)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$m <- factor(sample(1:4, n, TRUE), ordered = TRUE)
d$g <- factor(rep(1:30, each = 10))
eta <- 0.5 * d$x
p1 <- plogis(-0.7 - eta); p2 <- plogis(0.8 - eta)
P <- cbind(p1, p2 - p1, 1 - p2)
d$yo <- apply(P, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))

cat("\n### the stored X of yo ~ mo(m) + cs(m)\n")
fr <- frm(bf(yo ~ mo(m) + cs(m)), family = sratio(), data = d,
          dry_run = "frame")
lp <- fr$linpreds[["yo.mu"]]
X <- as.matrix(lp[["X"]])
cat("X columns: ", paste(colnames(X), collapse = ","), "\n", sep = "")
cat("X column sums of squares: ",
    paste(sprintf("%.6g", colSums(X^2)), collapse = ","), "\n", sep = "")
cat("cs labels: ",
    paste(vapply(lp[["cs"]], `[[`, "", "label"), collapse = ","),
    "\n", sep = "")
rnk <- function(A) {
  s <- sqrt(colSums(A^2)); s[!(s > 0)] <- 1
  qr(sweep(A, 2L, s, "/"))$rank
}
Z <- do.call(cbind, lapply(lp[["cs"]], `[[`, "vals"))
M <- cbind(1, X)
cat("rank(cbind(1, X)) = ", rnk(M), "  ncol = ", ncol(M), "\n", sep = "")
for (j in seq_len(ncol(Z))) {
  cat("  j=", j, " rank(cbind(M, Z[,1:j])) = ",
      rnk(cbind(M, Z[, seq_len(j), drop = FALSE])),
      "  needs >= ", rnk(M) + j, "\n", sep = "")
}
cat("qr() default tol (the rank test's tolerance) = ",
    formals(qr.default)$tol, "\n", sep = "")

cat("\n### the flat direction, measured as a profile\n")
f_mo <- suppressWarnings(frm(bf(yo ~ mo(m) + cs(m)), family = sratio(),
                             data = d))
f_cs <- frm(bf(yo ~ cs(m)), family = sratio(), data = d)
cat("logLik(mo(m) + cs(m)) = ", sprintf("%.9f", as.numeric(logLik(f_mo))),
    "  df = ", attr(logLik(f_mo), "df"), "\n", sep = "")
cat("logLik(cs(m))         = ", sprintf("%.9f", as.numeric(logLik(f_cs))),
    "  df = ", attr(logLik(f_cs), "df"), "\n", sep = "")
cat("difference = ", sprintf("%.6e",
      as.numeric(logLik(f_mo)) - as.numeric(logLik(f_cs))), "\n", sep = "")
V <- suppressWarnings(vcov(f_mo))
ev <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
cat("vcov eigenvalues of mo(m)+cs(m): ",
    paste(sprintf("%.4g", ev), collapse = " "), "\n", sep = "")
cat("NaN standard errors: ", sum(is.nan(sqrt(diag(V)))), " of ",
    nrow(V), "\n", sep = "")

cat("\n### the compat registry row for cs_pred x group:ordinal_cs\n")
tb <- as.data.frame(frm_compat(feature_a = "cs_pred()"))
sel <- tb[tb$feature_b == "group:ordinal_cs" |
            grepl("ordinal_cs", tb$feature_b), , drop = FALSE]
for (i in seq_len(nrow(sel))) {
  cat("  ", sel$feature_a[i], " x ", sel$feature_b[i], " = ",
      sel$status[i], "\n", sep = "")
  cat("    note: ", substr(sel$note[i], 1L, 400L), "\n", sep = "")
}

cat("\n### the importance refusal with cs()\n")
r <- tryCatch(frm(bf(yo ~ cs(m) + (1 | g)), family = sratio(), data = d,
                  importance = 50L),
              error = function(e) paste("REFUSED:", short(e)))
cat(if (is.character(r)) r else "FITTED (no refusal)", "\n")
cat("\nDONE ", TAG, "\n")
