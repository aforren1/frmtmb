# Review of lane quadseam, odds and ends, one build per process:
# 1. standard errors of three one-block fits (the Hessian is now read
#    from the rule taped at the estimate), and on the lane whether the
#    boundary probe's re-taped objective at the estimate equals the
#    fit's objective;
# 2. skew normal with (1 | g), default start, seeds 1..10 of the
#    acceptance stream of dev/quadseam-escorder.R: the optimizer vector
#    in hex (bitwise no-op check of the Laplace-side escape), the
#    escape record, logLik and seconds.
#   Rscript quadseam-review-misc.R <base|lane>
a <- commandArgs(TRUE)
arm <- a[1]
base <- "C:/Users/adf44/source/r/rellib-r20"
lane <- "C:/Users/adf44/source/r/wt-quadseam-lib"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (arm == "lane") c(lane, base, ul) else c(base, ul))
suppressMessages(library(frmtmb))
cat("lib:", arm, dirname(find.package("frmtmb")), "\n")
src <- readLines("quadseam-review-acc.R")
eval(parse(text = src[grep("^gen <- function", src):
                        (grep("^# exact marginal", src) - 1L)]))
for (cs in list(c("bern500", 2), c("bern2sd3", 4), c("pois200", 3),
                c("bern2sd3", 2))) {
  g <- gen(cs[1], as.integer(cs[2]))
  f <- suppressMessages(frm(fo[[g$fam]], family = frm_fam[[g$fam]],
                            data = g$d, quadrature = TRUE))
  se <- sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed))
  pr <- if (arm == "lane") {
    fn <- frmtmb:::se_probe_fn(f)
    sprintf("probe-obj %.3e", fn(f$opt$par) - f$opt$objective)
  } else ""
  cat(sprintf("SE %s %s %s se %s %s\n", cs[1], cs[2], arm,
              paste(sprintf("%.8g", se), collapse = ","), pr))
}
# skew normal, the setting of dev/quadseam-escorder.R
esc <- readLines("../quadseam-escorder.R")
cat("escorder header:", head(esc, 3), sep = "\n")
for (seed in 1:10) {
  set.seed(seed)
  n <- 200
  xs <- -abs(rnorm(n)) * 3
  y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
  dd <- data.frame(y = y, xs = xs, g = factor(rep(1:20, length.out = n)))
  fo2 <- bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1)
  m3 <- mean((y - mean(y))^3) / stats::sd(y)^3
  bad <- list(betad = c(log(stats::sd(y)), 2 * sign(m3) + 0.5 * m3))
  for (st in c("default", "stall")) {
    w <- character()
    t0 <- proc.time()[[3]]
    f <- withCallingHandlers(suppressMessages(
      frm(fo2, family = skew_normal(), data = dd, quadrature = TRUE,
          start = if (st == "stall") bad)),
      warning = function(x) {
        w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
      })
    tt <- proc.time()[[3]] - t0
    se <- f$opt[["stationary_escape"]]
    qs <- f$opt[["quad_settle"]]
    cat(sprintf(paste0("SKEW %2d %-7s %s ll %.10f esc %s settle %s ",
                       "warn %s secs %.2f par %s\n"), seed, st, arm,
                as.numeric(logLik(f)),
                if (is.null(se)) "-" else paste(format(se, digits = 4),
                                                collapse = "/"),
                if (is.null(qs)) "-" else paste(format(qs, digits = 3),
                                                collapse = "/"),
                if (length(w)) paste(w, collapse = ",") else "-",
                tt, paste(sprintf("%a", f$opt$par), collapse = ",")))
  }
}
cat("END\n")
