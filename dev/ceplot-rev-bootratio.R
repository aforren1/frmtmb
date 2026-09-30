# Reviewer check (lane ceplot): the bootstrap band against the Wald band
# at more refits and other seeds than the worker's 60 at seed 5, for the
# crossed unseen g:h (data seed 49) and the mm(by = ) new members (data
# seed 45). Refits 200 at bootstrap seeds 5, 6 and 7.
#   Rscript dev/ceplot-rev-bootratio.R > dev/ceplot-rev-log/bootratio.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-shapes.R")
f4 <- function(v) paste(sprintf("%.3f", v), collapse = " ")
ce <- function(o, cond, ...) {
  suppressWarnings(conditional_effects(o, "x", resolution = 3,
                                       re_formula = NULL, conditions = cond,
                                       ...))$x
}
w <- function(d) d$upper__ - d$lower__
fa <- shape_fits()$A$fit
set.seed(45)
d <- data.frame(x = rnorm(300), g1 = factor(sample(1:10, 300, TRUE)),
                g2 = factor(sample(1:10, 300, TRUE)))
fl <- rep(c("a", "b"), each = 5)
d$f1 <- factor(fl[d$g1])
d$f2 <- factor(fl[d$g2])
u <- rnorm(10, 0, 1)
d$y <- rnorm(300, 1 + 0.5 * d$x + 0.5 * (u[d$g1] + u[d$g2]), 0.5)
fm <- frm(bf(y ~ x + (1 | mm(g1, g2, by = cbind(f1, f2)))),
          family = gaussian(), data = d)
cases <- list(
  list("crossed g=1,h=1 unseen", fa, list(g = "1", h = "1")),
  list("crossed g=1,h=2 observed", fa, list(g = "1", h = "2")),
  list("mm by f1=a,f2=b unset", fm, list(f1 = "a", f2 = "b")),
  list("mm by f1=f2=a unset", fm, list(f1 = "a", f2 = "a")),
  list("mm by g1=2 seen, f2=b new", fm, list(g1 = "2", f1 = "a", f2 = "b")))
for (cs in cases) {
  wd <- ce(cs[[2]], cs[[3]])
  cat(cs[[1]], "| wald width", f4(w(wd)), "\n")
  for (s in 5:7) {
    bo <- tryCatch(ce(cs[[2]], cs[[3]], band = "boot", boot = 200, seed = s),
                   error = function(e) conditionMessage(e))
    if (is.character(bo)) {
      cat("   seed", s, "ERROR", bo, "\n")
      next
    }
    cat("   seed", s, "boot200 width / wald", f4(w(bo) / w(wd)),
        "| est inside boot band:",
        all(bo$lower__ <= wd$estimate__ & wd$estimate__ <= bo$upper__), "\n")
  }
}
