# rel084 release review, consolidator's semantic finding (a): the
# separation early return of check_convergence() leaves conv_classes unset,
# now reached also by ald_cens_sep(). Can a gp "boundary"/"field holds"
# message (gp_basis_message(), which reads conv_classes) appear beside a
# separation warning on a fit at a nonzero code?
#   SEPA  asym_laplace(), every row cens right, mu ~ gp(x, k = 5) on a field
#         at the white limit of its basis (ls 0.7 on [0, 10], n 60)
#   SEPB  bernoulli separated by w, plus gp(x, k = 5) (the pre-existing
#         trigger, separation_check())
# Each over seeds 1..6; prints code, conv_classes, and every condition.
#   REVLIB=r23|r22 Rscript dev/reviews/rel084-review-sepgp.R
source("dev/reviews/rel084-review-common.R")
options(width = 170)
field <- function(seed, n = 60) {
  set.seed(seed)
  x <- sort(runif(n, 0, 10))
  K <- exp(-outer(x, x, "-")^2 / (2 * 0.7^2)) + diag(1e-8, n)
  list(x = x, f = as.vector(t(chol(K)) %*% rnorm(n)))
}
show <- function(tag, r) {
  f <- r$v
  if (is.null(f)) { cat(tag, "ERROR", cls_line(r), substr(r$said[[length(r$said)]]$msg, 1, 200), "\n"); return() }
  cc <- if (is.environment(f$cache)) f$cache$conv_classes else "no-cache"
  cat(sprintf("%s | code %s | conv_classes %s | %s\n", tag, f$opt$convergence,
              if (is.null(cc)) "UNSET" else paste(cc, collapse = ","),
              paste(vapply(r$said, function(s) paste0(s$k, ":", paste(intersect(s$cls,
                c("frmtmb_separation", "frmtmb_gp_basis", "frmtmb_boundary_fit",
                  "frmtmb_limit_fit", "frmtmb_ald_local", "frmtmb_ald_no_se")),
                collapse = "/"), if (!length(intersect(s$cls, c("frmtmb_separation",
                "frmtmb_gp_basis", "frmtmb_boundary_fit", "frmtmb_limit_fit",
                "frmtmb_ald_local", "frmtmb_ald_no_se")))) s$cls[[1]] else ""), ""),
                collapse = " ")))
  for (s in r$said) if ("frmtmb_gp_basis" %in% s$cls)
    cat("    GP", s$k, ":", substr(gsub("\n", " ", s$msg), 1, 140), "...",
        if (grepl("hold", s$msg)) "[says the field holds]" else "", "\n")
}
for (s in 1:6) {
  fl <- field(9600 + s)
  d <- data.frame(x = fl$x, y = fl$f + 0.3 * rnorm(60), cc = "right")
  show(sprintf("SEPA seed %d", s),
       conds_of(frm(y | cens(cc) ~ gp(x, k = 5), family = asym_laplace(),
                    data = d)))
}
for (s in 1:6) {
  fl <- field(9700 + s)
  w <- rnorm(60)
  d <- data.frame(x = fl$x, w = w, y = as.integer(w > 0))
  show(sprintf("SEPB seed %d", s),
       conds_of(frm(y ~ w + gp(x, k = 5), family = bernoulli(), data = d)))
}
cat("DONE\n")
