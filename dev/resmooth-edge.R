# Lane wt-resmooth. The refusals the fix newly reaches: a t2() with an
# re margin at an unseen level, a genuinely PARTIAL re_formula beside a
# factor smooth, and simulate(newdata = ) at an unseen level.
#   RESMOOTH_LIB=base Rscript dev/resmooth-edge.R > dev/resmooth-edge-before.txt
#   Rscript dev/resmooth-edge.R > dev/resmooth-edge-after.txt
base <- identical(Sys.getenv("RESMOOTH_LIB"), "base")
.libPaths(c(if (!base) "C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb from:", find.package("frmtmb"), "\n")
msg <- function(e) substr(conditionMessage(e), 1, 220)
show <- function(lab, expr) {
  v <- tryCatch({ force(expr); "OK" }, error = function(e)
    paste("ERROR:", msg(e)))
  cat(sprintf("%-34s %s\n", lab, v))
}

set.seed(5)
n <- 300
d <- data.frame(x = runif(n), f = factor(rep(c("a", "b", "c"),
                                            length.out = n)),
                h = factor(rep(1:4, length.out = n)),
                g = factor(rep(1:10, length.out = n)))
d$y <- sin(2 * pi * d$x) + rnorm(10, 0, 0.5)[d$g] +
  rnorm(4, 0, 0.4)[d$h] + rnorm(3, 0, 0.3)[d$f] + rnorm(n, 0, 0.3)
nd <- d[1:4, ]
nd$g <- factor("99", levels = c(levels(d$g), "99"))

f2 <- suppressWarnings(frm(bf(y ~ t2(x, g, bs = c("cr", "re"))), data = d))
show("t2 new level", frm_linpred(f2, newdata = nd))
show("t2 new level, anl = TRUE",
     frm_linpred(f2, newdata = nd, allow_new_levels = TRUE))
show("t2 new level, re_formula = NA",
     frm_linpred(f2, newdata = nd, re_formula = NA))

# a GENUINELY partial formula: two bar terms, one named
fp <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5) + (1 | h) +
                               (1 | f)), data = d))
show("partial ~(1|h) beside an fs term",
     frm_linpred(fp, re_formula = ~ (1 | h)))
show("both bars named = NULL",
     frm_linpred(fp, re_formula = ~ (1 | h) + (1 | f)))
show("~(1|nosuch)", frm_linpred(fp, re_formula = ~ (1 | nosuch)))
show("simulate partial", simulate(fp, nsim = 2, re_formula = ~ (1 | h)))

ffs <- suppressWarnings(frm(bf(y ~ s(x, g, bs = "fs", k = 5)), data = d))
show("simulate(newdata, NA) new fs level",
     simulate(ffs, nsim = 2, newdata = nd, re_formula = NA))
show("simulate(newdata, NA) new fs level, anl",
     simulate(ffs, nsim = 2, newdata = nd, re_formula = NA,
              allow_new_levels = TRUE))
show("simulate(newdata, NA) no g column",
     simulate(ffs, nsim = 2, newdata = nd[, c("y", "x")], re_formula = NA))
fgg <- suppressWarnings(frm(bf(y ~ s(x) + (1 | g)), data = d))
show("simulate(newdata, NA) (1|g) new level",
     simulate(fgg, nsim = 2, newdata = nd, re_formula = NA))
show("simulate(newdata, NA) (1|g) no g column",
     simulate(fgg, nsim = 2, newdata = nd[, c("y", "x")], re_formula = NA))

# frm_bootstrap keeps its whole-model default: the smooths are redrawn
set.seed(1)
bp <- frmtmb:::sim_re_plan(ffs, NA, smooths = TRUE)
cat("bootstrap plan blocks:", paste(bp$blocks, collapse = ","),
    "| every:", isTRUE(bp$every), "\n")
sp <- frmtmb:::sim_re_plan(ffs, NA)
cat("simulate plan blocks:", paste(sp$blocks, collapse = ","),
    "| every:", isTRUE(sp$every), "\n")

# The guard with its positive condition ABSENT. `every = TRUE` is what
# the plan would report if the fs check were dropped now that no smooth
# block is redrawn, and then simulate(newdata = ) would pass
# allow_new_levels = TRUE to the design builder and answer an unseen
# level as the population curve while claiming it was a fresh level.
rspec <- ffs$spec$responses[[1L]]
show("design with the guard absent (nl_ok = TRUE)",
     frmtmb:::sim_newdata_design(ffs, rspec, nd, TRUE))
show("design with the guard present (nl_ok = FALSE)",
     frmtmb:::sim_newdata_design(ffs, rspec, nd, FALSE))
