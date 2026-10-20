# Lane illse2, fix round 1, finding 1: the REML log-likelihood of frmtmb
# and lme4 against the exact restricted log-likelihood (256 bits) at
# each one's own estimates, and the REML maximum, on
# test-se-reparam-reml.R's data (seed 1), degrees 3 to 7.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-fx-remlll.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-p3reml.R")
eval(parse(text = src[grep("^mp_logdet <- function", src):
                        (grep("^for [(]seed in seeds[)]", src) - 1L)]))
for (deg in 3:7) {
  set.seed(1)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.5)[d$g] +
    stats::rnorm(200, 0, 0.3)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                                "+ (1 | g)"))
  X <- stats::model.matrix(stats::as.formula(paste("~", paste(xn,
                                     collapse = " + "))), d)
  rp <- reml_parts(X, d$y, d$g)
  cst <- -(nrow(X) - ncol(X)) / 2 * log(2 * pi)
  ex <- function(th) cst - num(rp$nreml(m(th)))
  f <- suppressWarnings(frm(fo, data = d, REML = TRUE))
  l <- suppressWarnings(lme4::lmer(fo, data = d, REML = TRUE))
  thf <- unname(f$opt$par)
  thl <- log(c(sigma(l), attr(lme4::VarCorr(l)$g, "stddev")))
  o <- stats::optim(thf, function(t) -ex(t), method = "BFGS",
                    control = list(reltol = 1e-14))
  cat(sprintf(paste0("deg %d %s | frm logLik - exact at its point %+.2e | ",
                     "lme4 logLik - exact at its point %+.2e | ",
                     "exact max - frm %.2e | exact max - lme4 %.2e\n"),
              deg, arm, as.numeric(logLik(f)) - ex(thf),
              as.numeric(logLik(l)) - ex(thl), -o$value - ex(thf),
              -o$value - ex(thl)))
}
