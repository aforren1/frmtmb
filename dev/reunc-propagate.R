# Lane wt-reunc: predict(propagate_error = ), the user's decision of
# 2026-09-23, which renames `param_uncertainty` AND widens what it
# switches off.
#
# The old argument held the PARAMETERS at their estimates and still
# drew the group effects, on the reasoning that an effect is not a
# parameter. The new one holds both, because the axis the argument
# names is whether the error in the estimates is propagated, and the
# group effects are estimated too. So there are now four corners, and
# this checks that each is the thing it claims to be:
#
#   re_formula   propagate_error   what the interval carries
#   NULL         TRUE              noise + parameters + this group
#   NULL         FALSE             noise only, at this group's effect
#   NA           TRUE              noise + parameters, average group
#   NA           FALSE             noise only, average group
#
# Widths must order NA/FALSE <= NULL/FALSE and NULL/FALSE < NULL/TRUE,
# and the two FALSE corners must differ in CENTRE, since that is the
# thing re_formula cannot express.
#
#   Rscript dev/reunc-propagate.R <lib>
lib <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(lib)) lib <- "C:/Users/adf44/source/r/reunc2-lib"
.libPaths(unique(c(lib, "C:/Users/adf44/source/r/rellib-r3",
                   "C:/Users/adf44/AppData/Local/R/win-library/4.6")))
suppressMessages(library(frmtmb))
cat("frmtmb from", find.package("frmtmb"), "\n")

set.seed(4)
G <- 8; m <- 5
d <- data.frame(g = factor(rep(seq_len(G), each = m)), x = rnorm(G * m))
u <- rnorm(G, 0, 1.2)
d$y <- 1 + 0.5 * d$x + u[d$g] + rnorm(G * m, 0, 0.8)
fit <- frm(bf(y ~ x + (1 | g)), data = d)
# a level whose effect is FAR from the population, so that holding it
# at its estimate is visibly different from dropping it. Read from the
# fit's own b rather than through ranef()'s shape, which is not what
# is under test here.
bhat <- fit$estimates[["b"]]
far <- which.max(abs(bhat))
nd <- data.frame(x = 0, g = factor(levels(d$g)[far], levels = levels(d$g)))
cat("predicting level ", levels(d$g)[far], ", whose effect is ",
    sprintf("%.3f", bhat[far]), "\n", sep = "")

corner <- function(rf, pe) {
  set.seed(11)
  p <- predict(fit, newdata = nd, ndraws = 4000, re_formula = rf,
               propagate_error = pe)
  # unname(): p's own column names would otherwise be pasted onto
  # these, giving est.Estimate and width.Q97.5
  c(est = unname(p[1, "Estimate"]), width = unname(p[1, 4] - p[1, 3]))
}
g <- list(c("NULL", "TRUE"), c("NULL", "FALSE"),
          c("NA", "TRUE"), c("NA", "FALSE"))
res <- rbind(corner(NULL, TRUE), corner(NULL, FALSE),
             corner(NA, TRUE), corner(NA, FALSE))
rownames(res) <- vapply(g, function(z) paste(z, collapse = " / "), "")
cat("\nfour corners, 4000 draws, one seed\n")
print(round(res, 4))

# The two FALSE corners are the load-bearing comparison and four
# printed decimals cannot carry it: "identical" is a word this lane's
# M1 section spends pages distinguishing from "equal to a few ulps".
# Both directions are reported at full precision and in ulps.
ulpd <- function(a, b) abs(a - b) / (.Machine$double.eps * abs(a))
cat("\nthe two FALSE corners, at full precision\n")
cat("  width NULL/FALSE ", sprintf("%.17g", res[2, "width"]), "\n",
    "  width NA/FALSE   ", sprintf("%.17g", res[4, "width"]), "\n",
    "  identical(): ", identical(res[2, "width"], res[4, "width"]),
    "; apart by ", sprintf("%.2f", ulpd(res[2, "width"], res[4, "width"])),
    " ulp\n", sep = "")
cat("  centre difference         ",
    sprintf("%.17g", res[2, "est"] - res[4, "est"]), "\n",
    "  the group's fitted effect ", sprintf("%.17g", bhat[far]), "\n",
    "  apart by ",
    sprintf("%.2f", ulpd(res[2, "est"] - res[4, "est"], bhat[far])),
    " ulp\n", sep = "")

cat("\nthe orderings the argument claims\n")
cat("  NULL/FALSE narrower than NULL/TRUE: ",
    res[2, "width"] < res[1, "width"], "\n", sep = "")
cat("  NA/FALSE narrower than NA/TRUE:     ",
    res[4, "width"] < res[3, "width"], "\n", sep = "")
cat("  the two FALSE corners differ in CENTRE (what re_formula does,\n")
cat("  and propagate_error cannot): ",
    abs(res[2, "est"] - res[4, "est"]) > 0.5, " by ",
    sprintf("%.4f", abs(res[2, "est"] - res[4, "est"])), "\n", sep = "")
cat("  the two NULL corners agree in CENTRE (what propagate_error\n")
cat("  does, and re_formula cannot): ",
    abs(res[1, "est"] - res[2, "est"]) < 0.2, "\n", sep = "")

# The claim that FALSE now carries the observation noise ALONE. With
# sigma known that interval would be exactly 2 * 1.96 * sigma wide;
# here sigma is estimated, so the check is against the fit's own sigma
# rather than against a number typed in.
sig <- as.numeric(summary(fit)$spec_pars["sigma", "Estimate"])
cat("\nnoise-only check: NULL/FALSE width ",
    sprintf("%.4f", res[2, "width"]), " against 2 * 1.96 * sigma = ",
    sprintf("%.4f", 2 * 1.96 * sig), ", ratio ",
    sprintf("%.4f", res[2, "width"] / (2 * 1.96 * sig)), "\n", sep = "")

cat("\nthe old name is gone, not deprecated\n")
e <- tryCatch({
  predict(fit, newdata = nd, ndraws = 10, param_uncertainty = FALSE)
  "ACCEPTED"
}, error = function(e) conditionMessage(e))
cat("  predict(param_uncertainty = FALSE): ", substr(e, 1, 70), "\n",
    sep = "")
cat("  formals name it: ",
    "propagate_error" %in% names(formals(frmtmb:::predict.frmtmb_fit)),
    "; old name present: ",
    "param_uncertainty" %in% names(formals(frmtmb:::predict.frmtmb_fit)),
    "\n", sep = "")
