# Reviewer, punch round 2: what the mutant p2_mmsplit (two unset members
# drawn independently, weights not added) changes, on the lane's own
# mm(g1, g2) shape (data seed 43), draws (200 hand-built, seed 1) and
# boot (60 refits, seed 5), nothing set, before and after the mutation.
.libPaths(c("C:/Users/adf44/source/r/wt-postfit2-lib", "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
f3 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
set.seed(43)
d3 <- data.frame(x = rnorm(200), g1 = factor(sample(1:10, 200, TRUE)), g2 = factor(sample(1:10, 200, TRUE)))
u <- rnorm(10, 0, 1)
d3$y <- rnorm(200, 1 + 0.5 * d3$x + 0.5 * (u[d3$g1] + u[d3$g2]), 0.5)
f <- frm(bf(y ~ x + (1 | mm(g1, g2))), family = gaussian(), data = d3)
tpl <- f$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) f$estimates[[cp]]))
set.seed(1)
M <- matrix(rep(est, each = 200) + rnorm(200 * length(est), 0, 0.05), 200, dimnames = list(NULL, frmtmb::brms_par_labels(f)))
ds <- structure(list(stanfit = NULL, draws = cbind(frmtmb.sample:::draws_to_natural(M, f), lp__ = 0), fit = f), class = "frmtmb_draws")
go <- function(lab) {
  w <- suppressMessages(conditional_effects(f, "x", resolution = 3, re_formula = NULL))$x
  b <- suppressMessages(conditional_effects(f, "x", resolution = 3, re_formula = NULL, band = "boot", boot = 60, seed = 5))$x
  dd <- suppressMessages(conditional_effects(ds, "x", resolution = 3, re_formula = NULL, seed = 1))$x
  cat(lab, "widths: wald", f3(w$upper__ - w$lower__), "| boot", f3(b$upper__ - b$lower__), "| draws", f3(dd$upper__ - dd$lower__), "\n")
}
go("lane:   ")
fn <- get("ce_plan_part", asNamespace("frmtmb"))
tx <- deparse(fn)
tx <- sub("u <- unique(vals[isnew])", "u <- vals[isnew]", tx, fixed = TRUE)
tx <- sub("mt <- match(vals, u)", "mt <- cumsum(isnew)", tx, fixed = TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-postfit2/dev/postfit2-rev-p2-mmsplit-fix.R", local = TRUE)
nf <- eval(parse(text = tx)); environment(nf) <- asNamespace("frmtmb")
assignInNamespace("ce_plan_part", nf, "frmtmb")
go("mutant: ")
