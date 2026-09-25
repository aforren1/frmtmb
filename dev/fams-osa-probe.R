# residuals(type = "osa") on every core family with a zi, hu or zoi
# dpar, one small fit each.
#
#   FRMTMB_LIB=base Rscript dev/fams-osa-probe.R
#   FRMTMB_LIB=/opt/rlib/lane-fams Rscript dev/fams-osa-probe.R
#
# Output of both, convergence warnings dropped: dev/fams-osa-probe.txt.
lib <- Sys.getenv("FRMTMB_LIB", "/opt/rlib/lane-fams")
.libPaths(c(if (lib != "base") lib, "/opt/rlib/base", "/opt/rlib/deps",
            "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
set.seed(3)
n <- 120
x <- rnorm(n)
reg <- frmtmb:::family_registry
fams <- Filter(function(nm) {
  f <- tryCatch(reg[[nm]](), error = function(e) NULL)
  !is.null(f) && any(c("zi", "hu", "zoi") %in% f$dpars)
}, names(reg))
print(fams)
ys <- list(
  cnt = ifelse(runif(n) < 0.3, 0, rpois(n, 3) + 1),
  pos = ifelse(runif(n) < 0.3, 0, rgamma(n, 2, 1)),
  unit = ifelse(runif(n) < 0.2, 0,
                ifelse(runif(n) < 0.1, 1, rbeta(n, 2, 3))),
  line = ifelse(runif(n) < 0.3, 0, rnorm(n)))
pick <- function(nm) {
  if (grepl("beta", nm)) return("unit")
  if (grepl("gamma|lognormal", nm)) return("pos")
  if (grepl("laplace", nm)) return("line")
  "cnt"
}
for (nm in fams) {
  d <- data.frame(x = x, y = ys[[pick(nm)]], k = 10L)
  if (nm == "zero_inflated_beta") d$y[d$y == 1] <- 0.99
  form <- bf(y ~ x)
  if (nm == "zero_inflated_binomial") {
    form <- bf(y | trials(k) ~ x)
    d$y <- pmin(d$y, 10)
  }
  fit <- try(frm(form, family = reg[[nm]](), data = d), silent = TRUE)
  if (inherits(fit, "try-error")) {
    cat(nm, "FIT FAILED", conditionMessage(attr(fit, "condition")), "\n")
    next
  }
  r <- tryCatch({
    residuals(fit, type = "osa")
    "WORKS"
  }, error = function(e) conditionMessage(e))
  cat(sprintf("%-28s %s\n", nm, substr(r, 1, 110)))
}
