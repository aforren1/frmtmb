# Reviewer, item 1: confint()/vcov(full = TRUE) names on two more
# mixture shapes, hurdle x hurdle and a cs() component. Seed 20261006.
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261006)
n <- 500
d <- data.frame(x = rnorm(n), z = rnorm(n))
cls <- rbinom(n, 1, 0.4)
lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
d$yh <- ifelse(runif(n) < 0.15, 0L, d$y)
chk <- function(nm, fit) {
  rn <- rownames(confint(fit)); vc <- rownames(vcov(fit, full = TRUE))
  ok <- vapply(rn, function(p) {
    r <- tryCatch(suppressMessages(confint(fit, parm = p)), error = function(e) NULL)
    !is.null(r) && nrow(r) == 1L && identical(unname(r[1, ]), unname(confint(fit)[p, ]))
  }, NA)
  cat(sprintf("%-14s rows %d unique %s vcov-same %s parm-ok %d/%d\n  %s\n", nm, length(rn),
              !anyDuplicated(rn), identical(rn, vc), sum(ok), length(ok), paste(rn, collapse = ", ")))
}
f1 <- suppressWarnings(frm(bf(yh ~ x, hu1 ~ 1, hu2 ~ 1), data = d,
                           family = mixture(hurdle_cumulative(), hurdle_cumulative())))
chk("hurdle_hurdle", f1)
for (s in 1:10) {
  f2 <- tryCatch(suppressWarnings(frm(bf(y ~ x, mu2 ~ cs(z)), data = d,
                                      family = mixture(cumulative(), sratio()))),
                 error = function(e) NULL)
  if (!is.null(f2)) break
  d$z <- rnorm(n)
}
if (!is.null(f2)) chk("cs_component", f2) else cat("cs_component: no fit in 10 tries\n")
