# Reviewer re-check (lane ceplot, punch round 1): the B2 guard (the
# rename refuses a row with a grouping variable unset) with its
# condition ABSENT and PRESENT, on the fit (Wald, boot) and on draws.
#   Rscript dev/ceplot-rev-b2guard.R lane|base
# Data seeds 49, 41, 42, 43; boot seed 1 (20 refits); draws seed 1.
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/wt-ceplot-lib",
          "C:/Users/adf44/source/r/rellib-r4",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
source("C:/Users/adf44/source/r/frmtmb-wt-ceplot/dev/ceplot-rev-shapes.R")
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
f4 <- function(v) paste(sprintf("%.3f", v), collapse = " ")
three <- function(label, fit, cond) {
  ds <- hand(fit, 60)
  for (b in c("wald", "boot", "draws")) {
    r <- tryCatch(suppressWarnings(switch(b,
      wald = conditional_effects(fit, "x", resolution = 3, re_formula = NULL,
                                 conditions = cond),
      boot = conditional_effects(fit, "x", resolution = 3, re_formula = NULL,
                                 conditions = cond, band = "boot",
                                 boot = 20, seed = 1),
      draws = conditional_effects(ds, "x", resolution = 3, re_formula = NULL,
                                  conditions = cond, seed = 1)))$x,
      error = function(e) conditionMessage(e))
    cat(sprintf("%-44s %-5s %s\n", label, b, if (is.character(r))
      paste("REFUSED:", substr(r, 1, 70)) else
        paste("est", f4(r$estimate__), "| width", f4(r$upper__ - r$lower__))))
  }
}
fa <- shape_fits()$A$fit
three("A crossed, g=1 h=1 unseen", fa, list(g = "1", h = "1"))
three("A crossed, g=1 set, h unset", fa, list(g = "1"))
set.seed(41)
dt <- expand.grid(subj = factor(1:12), trt = factor(c("a", "b")), r = 1:5)
dt <- dt[!(dt$subj == "3" & dt$trt == "b"), ]
dt$x <- rnorm(nrow(dt))
dt$y <- rnorm(nrow(dt), 1 + 0.5 * dt$x + (dt$trt == "b") + rnorm(12)[dt$subj] +
                rnorm(24)[as.integer(interaction(dt$trt, dt$subj))], 0.5)
ft <- frm(bf(y ~ x + trt + (1 | subj) + (1 | trt:subj)), family = gaussian(),
          data = dt)
three("trt+(1|subj)+(1|trt:subj), trt=b subj=3", ft,
      list(trt = "b", subj = "3"))
three("  same, trt=b, subj unset", ft, list(trt = "b"))
three("  same, nothing set", ft, list())
set.seed(42)
de <- crossed_data(42)
de$f <- factor(ifelse(as.integer(de$g) <= 3, "a", "b"))
fe <- frm(bf(y ~ x + (1 | gr(g, by = f)) + (1 | f:h)), family = gaussian(),
          data = de)
three("gr(g, by=f) + (1|f:h), nothing set", fe, list())
three("  same, f=a, h unset", fe, list(f = "a"))
three("  same, f=a, g=1, h=2", fe, list(f = "a", g = "1", h = "2"))
