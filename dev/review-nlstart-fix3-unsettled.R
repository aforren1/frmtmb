# Review of lane nlstart, third short check: try to trigger
# `frmtmb_nl_multistart_unsettled` on real designs. Fits where the probe
# wins (review-nlstart-fix3-inv.R at 1e-8, fresh two_sat_rev and mm rows),
# under the default control, `check_se = "ignore"`, `se = FALSE`,
# `restarts = 0` and a small `grad_tol`. Prints every condition class.
# Usage: REV_ARM=after Rscript dev/review-nlstart-fix3-unsettled.R
source("dev/review-nlstart-common.R")
src <- readLines("dev/review-nlstart-fresh.R")
eval(parse(text = src[grep("^F <- list", src):(grep("^cat\\(\"ROWHEAD", src) - 1L)]))
cases <- list()
for (seed in 1:5) {
  cases[[length(cases) + 1L]] <- local({
    set.seed(seed * 37 + 11)
    n <- 300
    d <- data.frame(x = round(rnorm(n), 3))
    lat <- 1.5 * exp(0.6 * d$x)
    d$y <- ordered(cut(lat + rlogis(n), c(-Inf, 0.3, 1, 2, Inf), labels = 1:4))
    d$xs <- d$x * 1e-8
    list(lab = paste("exp cumulative 1e-8", seed), d = d, fam = cumulative(),
         fm = bf(y ~ b1 * exp(b2 * xs), b1 + b2 ~ 1, nl = TRUE))
  })
}
for (r in list(c("bernoulli", "two_sat_rev", 11), c("binomial", "two_sat_rev", 11),
               c("poisson", "two_sat_rev", 22), c("bernoulli", "mm", 55),
               c("binomial", "mm", 22), c("gamma", "logis3", 11))) {
  srcB <- if (r[2] == "logis3") {
    s2 <- readLines("dev/review-nlstart-fresh2.R")
    e <- new.env()
    eval(parse(text = s2[grep("^F <- list", s2):(grep("^cat\\(\"ROWHEAD", s2) - 1L)]), e)
    e
  } else environment()
  d <- srcB$gen(r[1], r[2], as.numeric(r[3]) * 1000 + match(r[2], names(srcB$B)))
  cases[[length(cases) + 1L]] <- list(lab = paste(r, collapse = " "), d = d,
                                      fam = F[[r[1]]], fm = srcB$form(r[1], r[2]))
}
ctls <- list(default = list(), check_se_ignore = list(control =
               frmtmb_control(check_se = "ignore")),
             se_false = list(se = FALSE),
             restarts0 = list(control = frmtmb_control(restarts = 0)),
             grad_tol_small = list(control = frmtmb_control(grad_tol = 1e-8)))
n_uns <- 0
for (cs in cases) for (cn in names(ctls)) {
  f <- tt(do.call(frm, c(list(cs$fm, family = cs$fam, data = cs$d), ctls[[cn]])))
  cls <- c(attr(f, "wc"), attr(f, "mc"))
  uns <- any(cls == "frmtmb_nl_multistart_unsettled")
  n_uns <- n_uns + uns
  cat(sprintf("%-32s %-16s ll %.4f conv %s | %s\n", cs$lab, cn, llv(f),
              if (iserr(f)) NA else f$opt$convergence,
              if (length(cls)) paste(cls, collapse = ",") else "-"))
  flush.console()
}
cat("unsettled raised:", n_uns, "of", length(cases) * length(ctls), "fits\n")
