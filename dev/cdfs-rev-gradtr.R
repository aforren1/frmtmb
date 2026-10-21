# Reviewer: gradients of TRUNCATED objectives (both bounds, and lower
# only) against Richardson differences, at the start and with the
# distributional intercept moved to +-4 on its link scale.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
set.seed(77)
n <- 200
x <- rnorm(n)
win <- function(y, lb, ub) { keep <- y >= lb & y <= ub; data.frame(x = x[keep], y = y[keep]) }
cases <- list(
  list("student", student(), win(1 + 0.5 * x + rt(n, 4), -1, 3), -1, 3),
  list("skew_normal", skew_normal(), win(1 + 0.5 * x + rnorm(n), 0, 3), 0, 3),
  list("Gamma", Gamma(link = "log"), win(rgamma(n, 3, 3 / exp(0.5 + 0.3 * x)), 0.5, 5), 0.5, 5),
  list("Beta", Beta(), win(rbeta(n, 2, 3), 0.1, 0.8), 0.1, 0.8),
  list("von_mises", von_mises(), win(rnorm(n, 0.3 * x, 0.6), -1.5, 1.5), -1.5, 1.5),
  list("negbinomial", negbinomial(), win(rnbinom(n, mu = exp(1 + 0.3 * x), size = 2), 1, 12), 1, 12),
  list("binomial", binomial(), cbind(win(rbinom(n, 12, plogis(0.3 * x)), 1, 10), tr = 12), 1, 10),
  list("hurdle_gamma", hurdle_gamma(), win(ifelse(runif(n) < 0.2, 0, rgamma(n, 3, 1.5)), 0, 5), NA, 5),
  list("exgaussian", exgaussian(), win(rnorm(n) + rexp(n), -1, 4), -1, 4),
  list("asym_laplace", asym_laplace(), win(rnorm(n), -1.5, 2), -1.5, 2)
)
for (cs in cases) {
  nm <- cs[[1]]; d <- cs[[3]]
  d$lbv <- cs[[4]]; d$ubv <- cs[[5]]
  tri <- nm == "binomial"
  dpn <- setdiff(frmtmb:::as_frmtmb_family(cs[[2]])$dpars, "mu")
  for (form in c("both", "lower")) {
    if (form == "lower" && is.na(cs[[4]])) next
    tt <- if (is.na(cs[[4]])) "trunc(ub = ubv)" else if (form == "both")
      "trunc(lb = lbv, ub = ubv)" else "trunc(lb = lbv)"
    if (form == "lower") d <- d  # rows already inside [lb, ub]
    lhs <- paste(c("y |", if (tri) "trials(tr) +", tt), collapse = " ")
    fo <- do.call(bf, c(list(as.formula(paste(lhs, "~ x"))),
                        lapply(dpn, function(dp) as.formula(paste(dp, "~ 1")))))
    ob <- tryCatch(frm(fo, family = cs[[2]], data = d, dry_run = "objective")$obj,
                   error = function(e) { cat(nm, form, "BUILD ERROR", conditionMessage(e), "\n"); NULL })
    if (is.null(ob)) next
    pts <- list(start = ob$par)
    jd <- which(names(ob$par) == "betad")
    for (j in seq_along(jd)) for (v in c(-4, 4)) {
      p <- ob$par; p[jd[j]] <- v
      pts[[sprintf("%s=%g", dpn[j], v)]] <- p
    }
    out <- vapply(names(pts), function(k) {
      p <- pts[[k]]
      f0 <- ob$fn(p)
      g <- as.numeric(ob$gr(p))
      fd <- numDeriv::grad(ob$fn, p, method = "Richardson")
      e <- max(abs(g - fd) / pmax(1, abs(fd)))
      sprintf("%s fn %.6g err %.1e%s", k, f0, e, if (is.finite(f0) && !(e < 1e-5)) " <--" else "")
    }, "")
    cat(sprintf("%-13s %-5s %s\n", nm, form, paste(out, collapse = " | ")))
  }
}
