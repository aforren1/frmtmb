# rel084 release review: a kept hole refit says its own boundary messages
# (nlband round 2), including gpselimit's gp message (it carries
# frmtmb_boundary_fit). Case 1: student() on gaussian band data, pole 1.3
# (nu runs to its limit). Case 2 (consolidator's semantic finding c): a
# nonlinear body with a pole and a gp() term at its basis's white limit,
# y = a / (w + b) + c, c ~ gp(x, k = 5) with ls 0.7 on [0, 10]. Each
# default fit is set beside the fit from the truth with nl_starts = FALSE;
# messages are listed in order with their counts.
#   REVLIB=r23|r22 Rscript dev/reviews/rel084-review-nlgp.R
source("dev/reviews/rel084-review-common.R")
options(width = 170)
ll <- function(r) if (is.null(r$v)) NA else as.numeric(logLik(r$v))
seq_of <- function(r) {
  if (!length(r$said)) return("-")
  paste(vapply(r$said, function(s) paste0(s$k, ":", s$cls[[1]],
    if ("frmtmb_boundary_fit" %in% s$cls) "(bnd)" else ""), ""), collapse = " ")
}
pole <- 1.3
form <- bf(y ~ a / (x + b), a + b ~ 1, nl = TRUE)
for (s in 1:3) {
  set.seed(970 + s)
  x <- runif(300, 0, 3)
  x <- x[abs(x - pole) > 0.25]
  d <- data.frame(x = x, y = 1 / (x - pole) + rnorm(length(x), 0, 0.2))
  r <- conds_of(frm(form, data = d, family = student()))
  q <- conds_of(frm(form, data = d, family = student(),
                    start = list(beta = c(1, -pole)),
                    control = frmtmb_control(nl_starts = FALSE)))
  cat(sprintf("STUDG seed %d default ll %.6f [%s]\n   truth  ll %.6f [%s]\n",
              s, ll(r), seq_of(r), ll(q), seq_of(q)))
}
formg <- bf(y ~ a / (w + b) + c, a + b ~ 1, c ~ gp(x, k = 5), nl = TRUE)
for (s in 1:3) {
  set.seed(9600 + s)
  n <- 60
  x <- sort(runif(n, 0, 10))
  K <- exp(-outer(x, x, "-")^2 / (2 * 0.7^2)) + diag(1e-8, n)
  f <- as.vector(t(chol(K)) %*% rnorm(n))
  w <- runif(400, 0, 3)
  w <- w[abs(w - pole) > 0.25][seq_len(n)]
  d <- data.frame(x = x, w = w, y = 1 / (w - pole) + f + 0.3 * rnorm(n))
  r <- conds_of(frm(formg, data = d))
  cat(sprintf("NLGP seed %d default ll %.6f b %s [%s]\n", s, ll(r),
              if (is.null(r$v)) NA else paste(format(
                r$v$opt$par[names(r$v$opt$par) == "beta"], digits = 4),
                collapse = " "),
              seq_of(r)))
  gpm <- msg_of(r, "frmtmb_gp_basis")
  if (!is.na(gpm)) cat("   GP MSG:", substr(gsub("\n", " ", gpm), 1, 260), "\n")
  q <- conds_of(frm(formg, data = d,
                    start = list(beta = c(1, -pole, if (!is.null(r$v))
                      r$v$opt$par[names(r$v$opt$par) == "beta"][-(1:2)]
                      else NULL)),
                    control = frmtmb_control(nl_starts = FALSE)))
  if (is.null(q$v)) q <- conds_of(frm(formg, data = d,
                    start = list(beta = c(1, -pole, 0)),
                    control = frmtmb_control(nl_starts = FALSE)))
  cat(sprintf("   truth  ll %.6f [%s]\n", ll(q), seq_of(q)))
}
cat("DONE\n")
