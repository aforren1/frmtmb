# Reviewer of lane fixes, claim 4c: confint() / vcov(full = TRUE) /
# vcov_cluster(full = TRUE) names over every ordinal family and
# threshold structure, cs(), grouped and multivariate thresholds.
#   Rscript dev/fixes-rev-confint.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("LIB", find.package("frmtmb"), "\n")
set.seed(31)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)),
                cl = factor(rep(1:30, 10)))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yg <- ifelse(d$g == "b", pmin(d$y, 4L), d$y)
d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
d$w <- 1L + (u + rnorm(n) > 0) + (u + rnorm(n) > 1)
fams <- c("cumulative", "sratio", "cratio", "acat")
cases <- list()
for (f in fams) for (th in c("flexible", "equidistant", "sum_to_zero")) {
  cases[[paste(f, th)]] <- list(bf(y ~ x), get(f)(threshold = th))
}
cases[["hurdle flexible"]] <- list(bf(yh ~ x), hurdle_cumulative())
cases[["hurdle equidistant"]] <- list(bf(yh ~ x),
                                      hurdle_cumulative(threshold = "equidistant"))
cases[["sratio cs"]] <- list(bf(y ~ x + cs(z)), sratio())
cases[["acat cs equidistant"]] <- list(bf(y ~ x + cs(z)),
                                       acat(threshold = "equidistant"))
cases[["cumulative grouped"]] <- list(bf(yg | thres(gr = g) ~ x), cumulative())
cases[["cumulative grouped equidistant"]] <-
  list(bf(yg | thres(gr = g) ~ x), cumulative(threshold = "equidistant"))
cases[["cumulative disc"]] <- list(bf(y ~ x, disc ~ 0 + z), cumulative())
cases[["mv cumulative + sratio"]] <- list(
  bf(y ~ x) + cumulative() + bf(w ~ x) + sratio(), NULL)
cases[["mv cumulative + cumulative equidistant"]] <- list(
  bf(y ~ x + cs(z)) + sratio() + bf(w ~ x + cs(z)) + sratio(), NULL)
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- tryCatch(suppressMessages(suppressWarnings(
    if (is.null(cs[[2]])) frm(cs[[1]], data = d) else
      frm(cs[[1]], family = cs[[2]], data = d))),
    error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("%-40s frm ERROR %s\n", nm, substr(conditionMessage(fit), 1,
                                                    120)))
    next
  }
  ci <- confint(fit)
  rn <- rownames(ci)
  V <- vcov(fit, full = TRUE)
  Vc <- tryCatch(vcov_cluster(fit, cluster = d$cl, full = TRUE),
                 error = function(e) conditionMessage(e))
  vc_ok <- if (is.character(Vc)) paste("vcov_cluster ERROR", substr(Vc, 1, 60))
    else identical(rownames(Vc), rn)
  # every row resolves by its own name to itself
  self <- vapply(seq_along(rn), function(i) {
    r <- tryCatch(suppressMessages(rownames(confint(fit, parm = rn[i]))),
                  error = function(e) "ERR")
    identical(r, rn[i])
  }, NA)
  tpl <- frmtmb:::outer_par_names(fit)
  old <- vapply(seq_along(tpl), function(i) {
    r <- tryCatch(suppressMessages(rownames(confint(fit, parm = tpl[i]))),
                  error = function(e) "ERR")
    identical(r, rn[i])
  }, NA)
  bfx <- grep("Intercept\\[", rn, value = TRUE)
  bfx <- bfx[!grepl("^log", bfx)]
  b_ok <- vapply(bfx, function(r) {
    identical(tryCatch(rownames(confint(fit, parm = paste0("b_", r))),
                       error = function(e) "ERR"), r)
  }, NA)
  cat(sprintf("%-40s unique %s  vcov %s  cluster %s  self %d/%d  old %d/%d  b_ %d/%d\n    %s\n",
              nm, !anyDuplicated(rn), identical(rownames(V), rn), vc_ok,
              sum(self), length(self), sum(old), length(old), sum(b_ok),
              length(b_ok),
              paste(rn, collapse = " | ")))
}
# the fixef() value of a threshold that IS the internal parameter
fit <- frm(bf(y ~ x), family = sratio(), data = d)
cat("sratio Intercept[3] confint est == fixef:",
    identical(unname(confint(fit, parm = "Intercept[3]")[, "est"]),
              unname(fixef(fit)["Intercept[3]", "Estimate"])), "\n")
# profile and uniroot by the new name
for (m in c("profile", "uniroot")) {
  r <- tryCatch(confint(fit, parm = "Intercept[2]", method = m),
                error = function(e) conditionMessage(e))
  cat(m, ":", if (is.character(r)) r else paste(rownames(r),
                                                format(r, digits = 6)), "\n")
}
fc <- frm(bf(y ~ x), family = cumulative(), data = d)
r <- tryCatch(confint(fc, parm = "log(Intercept[3] - Intercept[2])",
                      method = "profile"),
              error = function(e) conditionMessage(e))
cat("cumulative profile log-increment:", if (is.character(r)) r else
  paste(rownames(r), format(r, digits = 6)), "\n")
r <- tryCatch(confint(fc, parm = "Intercept[2]"),
              error = function(e) conditionMessage(e))
cat("cumulative parm = 'Intercept[2]' (a log increment internally):",
    if (is.character(r)) substr(r, 1, 300) else rownames(r), "\n")
