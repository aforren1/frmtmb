# Lane fixes, item 4b: brms 2.23.0's default_prior() rows for ordinal
# thresholds, and whether prior(class = "Intercept", coef = "k") is
# accepted, beside frmtmb's default_prior().
#   Rscript dev/fixes-prior-brms.R <lib>
lib <- commandArgs(TRUE)[1]
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb)})
cat("LIB", find.package("frmtmb"), "\n")
set.seed(62)
n <- 200
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)))
u <- stats::rlogis(n) + 0.8 * d$x
d$y <- 1L + (u > -1.2) + (u > -0.2) + (u > 0.8) + (u > 1.8)
d$yg <- ifelse(d$g == "b", pmin(d$y, 4L), d$y)
d$yh <- ifelse(runif(n) < 0.25, 0L, d$y)
show <- function(p) {
  p <- as.data.frame(p)
  keep <- p$class %in% c("Intercept", "delta")
  print(p[keep, intersect(c("prior", "class", "coef", "group", "resp",
                             "dpar", "lb", "ub", "source"), names(p))],
        row.names = FALSE)
}
cases <- list(
  flexible = list(y ~ x, "cumulative", list()),
  sratio_flexible = list(y ~ x, "sratio", list()),
  sum_to_zero = list(y ~ x, "cumulative", list(threshold = "sum_to_zero")),
  equidistant = list(y ~ x, "cumulative", list(threshold = "equidistant")),
  acat_equidistant = list(y ~ x, "acat", list(threshold = "equidistant")),
  grouped = list(yg | thres(gr = g) ~ x, "cumulative", list()),
  hurdle = list(yh ~ x, "hurdle_cumulative", list()))
for (nm in names(cases)) {
  cs <- cases[[nm]]
  cat("\n=====", nm, "\n")
  bfam <- do.call(get(cs[[2]], asNamespace("brms")), cs[[3]])
  ffam <- do.call(get(cs[[2]], asNamespace("frmtmb")), cs[[3]])
  cat("-- brms default_prior\n")
  show(brms::default_prior(brms::bf(cs[[1]]), data = d, family = bfam))
  cat("-- frmtmb default_prior\n")
  r <- tryCatch(default_prior(bf(cs[[1]]), data = d, family = ffam),
                error = function(e) e)
  if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n") else
    show(r)
  # brms's verdict on a per-threshold prior
  for (cf in c("1", "2")) {
    pr <- brms::prior_string("normal(0, 1)", class = "Intercept", coef = cf,
                             group = if (nm == "grouped") "a" else "")
    v <- tryCatch({
      brms::validate_prior(pr, brms::bf(cs[[1]]), data = d, family = bfam)
      "accepted"
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat("brms prior(Intercept, coef =", cf, "):", v, "\n")
    sc <- tryCatch({
      code <- brms::stancode(brms::bf(cs[[1]]), data = d, family = bfam,
                             prior = pr)
      grep("normal_lpdf\\(Intercept|normal_lpdf\\(first_Intercept",
           strsplit(code, "\n")[[1]], value = TRUE)
    }, error = function(e) paste("ERROR:", conditionMessage(e)))
    cat("   stancode line:", trimws(sc), "\n")
  }
}
