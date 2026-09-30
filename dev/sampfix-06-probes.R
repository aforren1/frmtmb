# Lane sampfix, script 06: three constructions around the laplace layout.
#
#   Rscript dev/sampfix-06-probes.R <lane|ref>
#
# (a) posterior_epred(newdata =) on laplace draws of y ~ x + (1 | g):
#     finite on the reference build; are the values right? Compared with
#     what the misaligned read would give, b[1] = log(sigma), b[2] =
#     theta_1, and with the population prediction.
# (b) a model whose only integrated component is mi() (no group-level
#     block): log_lik() and posterior_epred() on laplace draws.
# (c) laplace = TRUE on a REML fit: which labels land on which values.
# Data seed 1212 (as dev/arcovsample-rev-11-laplace.R), sampler seed 3.

arm <- commandArgs(trailingOnly = TRUE)[1L]
LANE <- "C:/Users/adf44/source/r/wt-sampfix-lib"
REF <- "C:/Users/adf44/source/r/rellib-r3"
USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (identical(arm, "lane")) c(LANE, REF, USER) else c(REF, USER))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
cat("ARM ", arm, "\n", sep = "")
show <- function(lab, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) e)
  if (inherits(r, "error")) {
    cat("  [", lab, "] ERROR: ", substr(conditionMessage(r), 1, 150), "\n",
        sep = "")
  } else {
    cat("  [", lab, "] ", paste(dim(r) %||% length(r), collapse = "x"),
        " non-finite ", sum(!is.finite(r)), "\n", sep = "")
  }
  invisible(r)
}
`%||%` <- function(a, b) if (is.null(a)) b else a

set.seed(1212L)
dd <- data.frame(g = factor(rep(1:6, each = 5L)), t = rep(1:5, 6L))
dd$x <- rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + rnorm(6L, 0, 0.5)[as.integer(dd$g)] +
  rnorm(nrow(dd), 0, 0.7)

cat("\n(a) newdata prediction on laplace draws\n")
fit <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd)
lap <- suppressWarnings(suppressMessages(
  frm_sample(fit, chains = 1, iter = 300, refresh = 0, seed = 3,
             laplace = TRUE)))
nd <- data.frame(x = c(-1, 1), g = factor(c(1, 2), levels = 1:6))
ep <- show("posterior_epred(newdata)", posterior_epred(lap, newdata = nd))
if (is.matrix(ep)) {
  m <- lap$draws
  pop <- cbind(m[, "b_Intercept"] - m[, "b_x"], m[, "b_Intercept"] + m[, "b_x"])
  cat("  max |epred - population|: ", format(max(abs(ep - pop))), "\n",
      sep = "")
  # which stored value each row's "group effect" was read from
  cand <- cbind(log_sigma = log(m[, "sigma"]), theta_1 = m[, "theta_1"],
                lp__ = m[, "lp__"])
  for (k in 1:2) {
    d <- apply(cand, 2L, function(v) max(abs(ep[, k] - pop[, k] - v)))
    cat("  row ", k, " (g = ", nd$g[k], "): epred - population equals ",
        names(which.min(d)), " to ", format(min(d)), "\n", sep = "")
  }
}
show("posterior_epred(newdata, re_formula = NA)",
     posterior_epred(lap, newdata = nd, re_formula = NA))

cat("\n(b) mi() as the only integrated component\n")
dm <- dd
dm$x[c(3, 11, 19)] <- NA
fm <- frm(bf(y ~ mi(x)) + bf(x | mi() ~ 1) + set_rescor(FALSE),
          family = gaussian(), data = dm)
cat("  template components: ",
    paste(names(fm$frame$par_template), collapse = " "), "\n", sep = "")
lm_ <- tryCatch(suppressWarnings(suppressMessages(
  frm_sample(fm, chains = 1, iter = 300, refresh = 0, seed = 3,
             laplace = TRUE))), error = function(e) e)
if (inherits(lm_, "error")) {
  cat("  frm_sample ERROR: ", conditionMessage(lm_), "\n", sep = "")
} else {
  cat("  columns: ", paste(colnames(lm_$draws), collapse = " "), "\n",
      sep = "")
  cat("  draws_is_laplace: ", frmtmb.sample:::draws_is_laplace(lm_), "\n",
      sep = "")
  show("log_lik", log_lik(lm_))
  show("posterior_epred resp y", posterior_epred(lm_, resp = "y"))
  show("posterior_epred resp x", posterior_epred(lm_, resp = "x"))
}

cat("\n(c) laplace = TRUE on a REML fit\n")
fr <- frm(bf(y ~ x + (1 | g)), family = gaussian(), data = dd, REML = TRUE)
cat("  random: ", paste(unique(names(fr$obj$env$par)[fr$obj$env$random]),
                        collapse = " "), "\n", sep = "")
for (pr in list(NULL, "flat")) {
cat(" prior = ", if (is.null(pr)) "default" else pr, "\n", sep = "")
lr <- tryCatch(suppressWarnings(suppressMessages(
  frm_sample(fr, chains = 1, iter = 300, refresh = 0, seed = 3,
             laplace = TRUE, prior = pr))), error = function(e) e)
if (inherits(lr, "error")) {
  cat("  frm_sample ERROR: ", substr(conditionMessage(lr), 1, 200), "\n",
      sep = "")
} else {
  cat("  columns: ", paste(colnames(lr$draws), collapse = " "), "\n",
      sep = "")
  cat("  column means: ", paste(format(colMeans(lr$draws), digits = 4),
                                collapse = " "), "\n", sep = "")
  cat("  stanfit parameter names: ",
      paste(unique(sub("[[].*", "", dimnames(as.array(lr$stanfit))[[3]])),
            collapse = " "), "\n", sep = "")
  cat("  ML fixef: ", paste(format(fixef(fr)[, 1], digits = 4),
                           collapse = " "), "\n", sep = "")
}
}
cat("\nDONE\n")
