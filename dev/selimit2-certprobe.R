# Lane selimit2 fix round: which zi/hu factor fit of test-se-limit.R's
# 4814 block loses its joint end under se_limit_joint_cert(), and why.
# Prints, per fit, the joint directions the search found before the
# certificate (se_limit_joint_dirs() at the stop), the row moves by
# outcome, and the certificate's answer.
#   Rscript dev/selimit2-certprobe.R
.libPaths(c("C:/Users/adf44/source/r/wt-selimit2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
n <- 300
for (b_end in c(FALSE, TRUE)) {
  set.seed(21)
  d <- data.frame(f = factor(rep(c("a", "b"), each = n / 2)))
  if (!b_end) d$x <- stats::rnorm(n)
  d$y <- ifelse(d$f == "a", 0L, stats::rpois(n, 4))
  for (fam in list(zero_inflated_poisson(), hurdle_poisson())) {
    dp <- if (identical(fam$family, "hurdle_poisson")) "hu" else "zi"
    fm <- stats::as.formula(paste(dp, "~ f"))
    ws <- character(0)
    fit <- withCallingHandlers(frm(bf(y ~ 1, fm), family = fam, data = d),
      warning = function(w) {
        ws <<- c(ws, substr(conditionMessage(w), 1, 120))
        invokeRestart("muffleWarning")
      }, message = function(m) invokeRestart("muffleMessage"))
    p <- fit$opt$par
    cat(sprintf("b_end %s %s | par %s | lost %s | warn %s\n", b_end, dp,
                paste(format(p, digits = 5), collapse = " "),
                paste(names(ns$sdr_of(fit)$se_lost), ns$sdr_of(fit)$se_lost,
                      collapse = ","), paste(ws, collapse = " || ")))
    om <- ns$outer_par_map(fit)
    u <- ns$se_limit_units(fit, om)
    H <- fit$obj$he(p)
    np <- length(p)
    cand <- ns$se_limit_candidates(fit)
    jd <- ns$se_limit_joint_dirs(fit, H, p, u, cand, om)
    fn <- function(q) fit$obj$fn(q)
    f0 <- fn(p)
    for (w in jd) {
      lp <- fit$frame$linpreds[[2]]
      pos <- ns$se_lp_outer_pos(fit, lp, om)
      a <- as.numeric(lp$X %*% w[pos])
      s <- ns$se_limit_joint_dir(fn, f0, p, w, H, attr(w, "z"))
      cat("  dir", format(w, digits = 4), "| s", s, "| a by group/y0:",
          paste(tapply(a, paste(d$f, d$y == 0), function(v)
            format(range(v), digits = 3)), collapse = " "), "\n")
      if (s != 0) {
        cat("  cert", ns$se_limit_joint_cert(fit, fn, f0, p, s * w, om), "\n")
      }
    }
  }
}
