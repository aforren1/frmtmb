# rev-gddm: is the recomputation noise that justifies the tolerance real,
# and how wide is the band the tolerance opens?
#
# Seed 202609, the sweep's own seed and its own design, so the number
# under test is the one dev/gddm-findings.md records (6.03e-15 on
# poly(cohn, 2), 391 ulp at the column maximum 0.06944).

lib <- Sys.getenv("GDDM_LIB", "C:/Users/adf44/source/r/rev-gddm-lib")
.libPaths(c(lib,
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm:", lib, " frmtmb.eam", format(packageVersion("frmtmb.eam")), "\n\n")

seed <- 202609
set.seed(seed)
n <- 480L
d <- gddm_simulate(n, mu = 2, bs = 2.5, ndt = 0.25,
                   control = gddm_control(t_max = 2))
d$subj <- factor(sample.int(6L, n, replace = TRUE))
d$cohn <- rep(c(0, 0.128, 0.256, 0.512), each = n / 4L)
d$coh <- factor(d$cohn)

cat("== 1. is poly() row-deterministic on identical inputs?\n")
p <- poly(d$cohn, 2)
gi <- match(d$cohn, sort(unique(d$cohn)))
first <- match(seq_len(4L), gi)
for (k in 1:2) {
  v <- p[, k]
  ref <- v[first][gi]
  dd <- abs(v - ref)
  cat(sprintf("  poly column %d: max |row - first row of its own value| = %.3e",
              k, max(dd)))
  cat(sprintf("   bitwise identical to first row on all %d rows: %s\n",
              n, identical(v, unname(ref))))
}
cat("  whole matrix, all rows vs their own value's first row: ",
    "max diff ", format(max(abs(p - p[first, ][gi, ])), digits = 3), "\n",
    sep = "")

cat("\n== 2. the same question the check asks: condition = coh\n")
cat("  (the sweep's poly design has one condition per coherence level)\n")
tolk <- vapply(1:2, function(k) 1e-8 * max(abs(p[, k])), 0)
cat(sprintf("  column max: %.5g %.5g   tolerance: %.3e %.3e\n",
            max(abs(p[, 1])), max(abs(p[, 2])), tolk[1], tolk[2]))

cat("\n== 3. does an EXACT comparison false-alarm on this design?\n")
d$cond <- gddm_conditions(d, cohn)
fit_try <- function(mu_rhs, data) {
  tryCatch({
    frm(bf(as.formula(paste("rt | vint(upper, cond) ~", mu_rhs)),
           bs ~ 1, ndt ~ 1, bias = 0.5),
        family = gddm(control = gddm_control(t_max = 2, dt = 0.05,
                                             ny = 51L)),
        data = data, dry_run = "frame")
    ""
  }, error = function(e) conditionMessage(e))
}
m <- fit_try("poly(cohn, 2)", d)
cat("  poly(cohn, 2) with cond = cohn, current tolerance: ",
    if (nzchar(m)) paste("REFUSED:", substr(m, 1, 90)) else "accepted", "\n",
    sep = "")

# What an exact comparison would say, computed on the frame the check
# itself reads, without changing the package: replay gd_varying_groups
# with a zero tolerance.
ne_exact <- function(v, gi, first) {
  m <- as.matrix(v)
  r <- m[first[gi], , drop = FALSE]
  bad <- apply(m != r, 1L, any)
  sort(unique(gi[which(bad)]))
}
cat("  exact comparison on the poly matrix, per condition: ",
    length(ne_exact(p, gi, first)), " conditions flagged\n", sep = "")

cat("\n== 4. ulp of the column maximum\n")
for (k in 1:2) {
  mx <- max(abs(p[, k]))
  cat(sprintf("  column %d max %.6g, 1 ulp = %.3e, 6.03e-15 = %.1f ulp\n",
              k, mx, .Machine$double.eps * mx, 6.03e-15 /
              (.Machine$double.eps * mx / 2)))
}

cat("\n== 5. other computed model-frame columns: are they row-deterministic?\n")
x <- d$cohn
probe <- list(
  "scale(x)"      = scale(x),
  "poly(x, 2)"    = poly(x, 2),
  "poly(x, 3)"    = poly(x, 3),
  "splines::bs(x, 4)" = splines::bs(x, 4),
  "splines::ns(x, 3)" = splines::ns(x, 3),
  "log1p(x)"      = log1p(x),
  "I(x/sd(x))"    = x / sd(x),
  "cut(x, 3)"     = as.integer(cut(x, 3)))
for (nm in names(probe)) {
  v <- as.matrix(probe[[nm]])
  r <- v[first[gi], , drop = FALSE]
  cat(sprintf("  %-20s max |row - first| = %.3e   bitwise equal: %s\n",
              nm, max(abs(v - r)), identical(as.vector(v), as.vector(r))))
}
