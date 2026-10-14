# Review of lane ci701: two siblings where the lane's scan saw the
# pull-back run. test-se-check.R RB4 ((1 + x | g2), 45 outer
# parameters, seeds 11 and 17) against lme4; test-mo-factor-re.R:102
# ((mo(x) | g), code 1, through se_boundary_stop()). Lane build, as
# installed and with se_edge_cor() stubbed, in one process.
#
#   Rscript dev/ci701-rev-sib.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
orig <- get("se_edge_cor", ns)
stub <- function(fit, np) logical(np)
cap <- function(expr) {
  out <- character(0)
  v <- withCallingHandlers(expr, warning = function(w) {
    out <<- c(out, paste("W:", substr(conditionMessage(w), 1, 90)))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    out <<- c(out, paste("M:", substr(conditionMessage(m), 1, 90)))
    invokeRestart("muffleMessage")
  })
  list(value = v, conds = out)
}
for (s in c(11L, 17L)) {
  set.seed(s)
  n <- 480
  d <- data.frame(x = rnorm(n), f = factor(sample(1:40, n, TRUE)),
                  g2 = factor(rep(1:6, length.out = n)))
  d$y <- 1 + 0.5 * d$x + rnorm(40, 0, 0.5)[d$f] + rnorm(n)
  m <- suppressMessages(suppressWarnings(
    lme4::lmer(y ~ x + f + (1 + x | g2), data = d, REML = FALSE)))
  ref <- sqrt(diag(as.matrix(stats::vcov(m))))
  for (a in c("lane", "noD2")) {
    assignInNamespace("se_edge_cor", if (a == "lane") orig else stub, ns)
    r <- cap({
      fit <- frm(y ~ x + f + (1 + x | g2), data = d)
      se <- sqrt(diag(vcov(fit)))
      list(fit = fit, se = se)
    })
    fit <- r$value$fit
    se <- r$value$se
    lost <- ns$sdr_of(fit)$se_lost
    stopifnot(length(se) == length(ref)); rr <- unname(se) / unname(ref) - 1; rn <- unname(se) / unname(ref)
    cat(sprintf(paste0("RES rb4 s=%d %s code=%d t=%.4g lost=%s ",
                       "maxrel_vs_lme4=%.3g medrel=%.3g\n"), s, a,
                fit$opt$convergence, fit$estimates$theta[3],
                paste(names(lost), lost, sep = ":", collapse = ","),
                max(abs(rr)), stats::median(abs(rr))))
    cat("   se[1:3]:", sprintf("%.6f", se[1:3]), " lme4:",
        sprintf("%.6f", ref[1:3]), "\n")
    for (x in r$conds) cat("   ", x, "\n")
  }
}
src <- readLines(file.path("C:/Users/adf44/source/r/frmtmb-wt-ci701",
                           "tests/testthat/test-mo-factor-re.R"))
i0 <- grep("^mo_fac_data <- function", src)
i1 <- i0 + which(src[(i0 + 1):length(src)] == "}")[1]
eval(parse(text = src[i0:i1]))
dd <- mo_fac_data()
for (a in c("lane", "noD2")) {
  assignInNamespace("se_edge_cor", if (a == "lane") orig else stub, ns)
  r <- cap(frm(bf(y ~ mo(x) + (mo(x) | g)), family = gaussian(), data = dd))
  fit <- r$value
  lost <- ns$sdr_of(fit)$se_lost
  se <- suppressWarnings(sqrt(diag(vcov(fit, full = TRUE))))
  cat(sprintf("RES mo %s code=%d t=%.4g lost=%s\n", a, fit$opt$convergence,
              fit$estimates$theta[3],
              paste(names(lost), lost, sep = ":", collapse = ",")))
  cat("   se:", sprintf("%.6f", se), "\n")
  for (x in r$conds) cat("   ", x, "\n")
}
assignInNamespace("se_edge_cor", orig, ns)
