# Reviewer, claim 5: cat() against brms 2.23.0. Seed 505.
# Log: dev/aterms2-rev-log-05-cat.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cond <- function(label, expr) {
  w <- character(0)
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
  cat(sprintf("%-40s %s\n", label,
              if (is.character(r) && length(r) == 1L) substr(r, 1, 160)
              else paste("ok:", class(r)[1])))
  for (x in unique(w)) cat("      warning:", substr(x, 1, 150), "\n")
  invisible(r)
}
set.seed(505)
n <- 150
d <- data.frame(x = rnorm(n))
d$o <- factor(cut(d$x + rnorm(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE),
              ordered = TRUE)
d$yi <- as.integer(d$o)
d$k <- 4
d$k5 <- 5
d$yg <- rnorm(n)

cat("brms family $ad containing 'cat':\n")
for (fm in c("cumulative", "sratio", "cratio", "acat", "categorical",
             "multinomial", "gaussian", "poisson", "binomial")) {
  ad <- tryCatch(brms:::brmsfamily(fm)$ad, error = function(e) "?")
  cat(sprintf("  %-12s %s\n", fm, paste(ad, collapse = ",")))
}
cat("\nbrms standata:\n")
b1 <- cond("brms cat(4)", brms::standata(brms::bf(yi | cat(4) ~ x), d,
                                         family = brms::cumulative()))
cat("   nthres", b1$nthres, "\n")
b2 <- cond("brms cat(5), 4 observed", brms::standata(
  brms::bf(yi | cat(5) ~ x), d, family = brms::cumulative()))
cat("   nthres", b2$nthres, "\n")
b3 <- cond("brms cat(k) column", brms::standata(
  brms::bf(yi | cat(k) ~ x), d, family = brms::cumulative()))
cat("   nthres", if (is.list(b3)) b3$nthres else NA, "\n")
b4 <- cond("brms cat(4) + thres(2)", brms::standata(
  brms::bf(yi | cat(4) + thres(2) ~ x), d, family = brms::cumulative()))
cat("   nthres", if (is.list(b4)) b4$nthres else NA, "\n")
cond("brms cat(3) on gaussian", brms::standata(
  brms::bf(yg | cat(3) ~ x), d, family = gaussian()))
cond("brms cat(4) on categorical", brms::standata(
  brms::bf(yi | cat(4) ~ x), d, family = brms::categorical()))

cat("\nfrmtmb:\n")
f1 <- cond("frm cat(4)", frm(yi | cat(4) ~ x, data = d,
                             family = cumulative()))
f0 <- cond("frm thres(3)", frm(yi | thres(3) ~ x, data = d,
                               family = cumulative()))
if (inherits(f1, "frmtmb_fit")) {
  set.seed(1)
  p <- f0$opt$par + rnorm(length(f0$opt$par), 0, 0.1)
  cat("   fn identical cat(4) vs thres(3):",
      identical(f1$obj$fn(p), f0$obj$fn(p)), "\n")
}
f2 <- cond("frm cat(5), 4 observed", frm(yi | cat(5) ~ x, data = d,
                                         family = cumulative()))
f2b <- cond("frm thres(4), 4 observed", frm(yi | thres(4) ~ x, data = d,
                                           family = cumulative()))
if (inherits(f2, "frmtmb_fit") && inherits(f2b, "frmtmb_fit")) {
  cat("   thresholds", length(f2$frame$par_template$tau_raw %||% 0),
      " logLik identical:", identical(logLik(f2), logLik(f2b)), "\n")
}
f3 <- cond("frm cat(k) column", frm(yi | cat(k) ~ x, data = d,
                                    family = cumulative()))
if (inherits(f3, "frmtmb_fit")) {
  cat("   fn identical to cat(4):", identical(f3$obj$fn(p), f1$obj$fn(p)),
      "\n")
}
cond("frm cat(4) + thres(2)", frm(yi | cat(4) + thres(2) ~ x, data = d,
                                  family = cumulative()))
cond("frm thres(2) + cat(4)", frm(yi | thres(2) + cat(4) ~ x, data = d,
                                  family = cumulative()))
cond("frm cat(3) on gaussian", frm(yg | cat(3) ~ x, data = d))
cond("frm cat(4) on categorical", frm(yi | cat(4) ~ x, data = d,
                                      family = categorical()))
cond("frm cat(4) on factor response", frm(o | cat(4) ~ x, data = d,
                                          family = cumulative()))
cond("frm cat(x = 4)", frm(yi | cat(x = 4) ~ x, data = d,
                           family = cumulative()))
cond("frm cat() empty", frm(yi | cat() ~ x, data = d,
                            family = cumulative()))
