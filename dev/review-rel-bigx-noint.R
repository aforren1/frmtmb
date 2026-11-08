# Release review 0.80.0: bigx's centered rank test (findings section 18)
# on designs with no "(Intercept)" column, where design_center() returns
# the design uncentered: an ordinal family (its thresholds replace the
# intercept) and `0 + f`. Product and power of columns shifted by `m`,
# against the same model on the centered columns.
#
#   Rscript dev/review-rel-bigx-noint.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(strsplit(lib, ";")[[1]], "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")),
    as.character(packageVersion("frmtmb")), "\n")
run <- function(expr) {
  cs <- character()
  v <- tryCatch(withCallingHandlers(expr,
    warning = function(w) {
      cs <<- c(cs, paste0("W:", class(w)[1L]))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cs <<- c(cs, paste0("M:", substr(conditionMessage(m), 1, 60)))
      invokeRestart("muffleMessage")
    }), error = function(e) {
      cs <<- c(cs, paste0("E:", substr(conditionMessage(e), 1, 80)))
      NULL
    })
  list(v = v, cs = unique(cs))
}
ll <- function(f) if (is.null(f)) NA else as.numeric(logLik(f))
for (seed in 1:3) for (m in c(1e4, 1e5, 1e6)) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n),
                  f = factor(sample(letters[1:3], n, TRUE)))
  eta <- 0.5 * d$x1 - 0.4 * d$x2 + 0.6 * d$x1 * d$x2 + 0.3 * d$x1^2
  d$yo <- cut(eta + rlogis(n), c(-Inf, -0.5, 0.5, 1.5, Inf),
              labels = 1:4, ordered_result = TRUE)
  d$yg <- c(a = 0, b = 1, c = -1)[as.character(d$f)] + eta + rnorm(n, 0, 0.5)
  d$w1 <- d$x1 + m
  d$w2 <- d$x2 + m
  cases <- list(
    cum_prod = list(bf(yo ~ w1 * w2), bf(yo ~ x1 * x2), cumulative()),
    cum_pow = list(bf(yo ~ w1 + I(w1^2)), bf(yo ~ x1 + I(x1^2)),
                   cumulative()),
    nof_prod = list(bf(yg ~ 0 + f + w1 * w2), bf(yg ~ 0 + f + x1 * x2),
                    gaussian()),
    icpt_prod = list(bf(yg ~ f + w1 * w2), bf(yg ~ f + x1 * x2),
                     gaussian()))
  for (nm in names(cases)) {
    cz <- cases[[nm]]
    a <- run(frm(cz[[1]], data = d, family = cz[[3]]))
    b <- run(frm(cz[[2]], data = d, family = cz[[3]]))
    kept <- if (!is.null(a$v)) length(a$v$estimates$beta) else NA
    cat(sprintf("seed %d m %.0e %-9s shifted ll %.4f (beta %s) centered %.4f (beta %d) short %.3f | %s\n",
                seed, m, nm, ll(a$v), kept, ll(b$v),
                length(b$v$estimates$beta), ll(b$v) - ll(a$v),
                if (length(a$cs)) paste(a$cs, collapse = " ; ") else "-"))
  }
}
cat("NOINT DONE\n")
