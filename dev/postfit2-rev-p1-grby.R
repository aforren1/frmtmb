# Reviewer, punch round 1: gr(g, by = f), f constant within g.
#   Rscript dev/postfit2-rev-p1-grby.R <base|lane>
args <- commandArgs(TRUE); arm <- if (length(args)) args[1] else "lane"
libs <- c("C:/Users/adf44/source/r/wt-postfit2-lib", "C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "base") libs <- libs[-1]
.libPaths(libs)
suppressMessages({library(frmtmb); library(frmtmb.sample)})
say <- function(...) cat(sprintf(...), "\n", sep = "")
f3 <- function(v) paste(sprintf("%.4f", v), collapse = " ")
say("ARM %s", arm)
set.seed(44)
d4 <- data.frame(x = rnorm(240), g = factor(rep(1:12, 20)))
d4$f <- factor(ifelse(as.integer(d4$g) <= 6, "a", "b"))
sdg <- ifelse(1:12 <= 6, 0.3, 1.5)
d4$y <- rnorm(240, 1 + 0.5 * d4$x + rnorm(12, 0, sdg)[d4$g], 0.5)
f4 <- frm(bf(y ~ x + f + (1 | gr(g, by = f))), family = gaussian(), data = d4)
bk <- f4$frame$re_blocks
say("blocks %d; group_name %s; levels %s", length(bk),
    paste(vapply(bk, function(b) b$group_name %||% "", ""), collapse = ","),
    paste(head(bk[[1]]$levels, 4), collapse = ","))
print(VarCorr(f4))
tpl <- f4$frame$par_template
est <- unlist(lapply(names(tpl), function(cp) f4$estimates[[cp]]))
set.seed(1)
M <- matrix(rep(est, each = 200) + rnorm(200 * length(est), 0, 0.05), 200,
            dimnames = list(NULL, frmtmb::brms_par_labels(f4)))
ds <- structure(list(stanfit = NULL, draws = cbind(
  frmtmb.sample:::draws_to_natural(M, f4), lp__ = 0), fit = f4), class = "frmtmb_draws")
run <- function(label, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) e)
  if (inherits(r, "error")) return(say("%s: ERROR %s", label, substr(conditionMessage(r), 1, 200)))
  say("%s: est %s | width %s", label, f3(r[[1]]$estimate__), f3(r[[1]]$upper__ - r[[1]]$lower__))
}
for (fl in c("a", "b")) {
  cnd <- list(f = fl)
  run(paste("f =", fl, "g unset | wald"), conditional_effects(f4, "x", resolution = 3, re_formula = NULL, conditions = cnd))
  run(paste("f =", fl, "g unset | boot40"), conditional_effects(f4, "x", resolution = 3, re_formula = NULL, conditions = cnd, band = "boot", boot = 40, seed = 5))
  run(paste("f =", fl, "g unset | draws"), conditional_effects(ds, "x", resolution = 3, re_formula = NULL, conditions = cnd, seed = 1))
}
for (cnd in list(list(f = "a", g = "3"), list(f = "b", g = "9"))) {
  lab <- sprintf("f = %s g = %s (observed)", cnd$f, cnd$g)
  run(paste(lab, "| wald"), conditional_effects(f4, "x", resolution = 3, re_formula = NULL, conditions = cnd))
  run(paste(lab, "| boot40"), conditional_effects(f4, "x", resolution = 3, re_formula = NULL, conditions = cnd, band = "boot", boot = 40, seed = 5))
  run(paste(lab, "| draws"), conditional_effects(ds, "x", resolution = 3, re_formula = NULL, conditions = cnd, seed = 1))
}
