# Re-check of lane illse2, finding 2: the escape rule optctrl_stops_earlier() (from dev/illse2-rev-esc2.R).
# controls" the right condition for sd_up_escape()? Degree 7 seed 4 of
# dev/illse2-p2.R with controls that ask for MORE work than the default.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-rev-esc2.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
deg <- 7L
set.seed(4)
n <- 300
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(1:20, length.out = n)))
d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
xn <- paste0("x", seq_len(deg))
for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                              "+ (1 | g)"))
# exact profile maximum from dev/illse2-findings.md section 2
llmax <- -61.7559980476
ctls <- list(
  more_iter_5000 = frmtmb_control(optCtrl = list(iter.max = 5000, eval.max = 5000)),
  looser_rel_tol = frmtmb_control(optCtrl = list(rel.tol = 1e-6)),
  lower_iter_50 = frmtmb_control(optCtrl = list(iter.max = 50, eval.max = 1000)),
  optim_BFGS = frmtmb_control(optimizer = "optim"),
  zero_iter = frmtmb_control(optCtrl = list(iter.max = 0)),
  default = frmtmb_control(),
  same_as_default = frmtmb_control(optCtrl = list(iter.max = 1000,
                                                  eval.max = 1000)),
  more_iter = frmtmb_control(optCtrl = list(iter.max = 5000,
                                            eval.max = 5000)),
  tighter_tol = frmtmb_control(optCtrl = list(iter.max = 1000,
                                              eval.max = 1000,
                                              rel.tol = 1e-12)),
  restarts3 = frmtmb_control(restarts = 3),
  grad_tol_1e4 = frmtmb_control(grad_tol = 1e-4)
)
for (nm in names(ctls)) {
  t0 <- proc.time()[[3]]
  r <- cap(frm(fo, data = d, control = ctls[[nm]]))
  tt <- proc.time()[[3]] - t0
  f <- r$v
  if (is_err(f)) { cat(nm, "ERROR", f, "\n"); next }
  p <- f$opt$par
  se <- suppressWarnings(sqrt(diag(vcov(f))))
  cat(sprintf("%-16s code %d | short %.3e | sd %.4g | escape %s | warn %s | lost %s | SE finite %d of %d | %.1f s\n",
              nm, f$opt$convergence, llmax - as.numeric(logLik(f)),
              exp(p[names(p) == "theta"]),
              if (is.null(f$opt$sd_escape)) "-" else
                format(f$opt$sd_escape[[1]], digits = 3),
              wtag(r$w), lost_of(f), sum(is.finite(se)), length(se), tt))
}
