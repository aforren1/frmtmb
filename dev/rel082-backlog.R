# The backlog after the union merge of the twelve lanes: their sections
# "Filed in the round of 2026-11-11" already sit as one section with one
# "### From lane <lane>" subsection each, above "Filed at the 0.81.0
# release", as the convention of the earlier release sections has it.
# This script renames that section, restores the blank line that lane
# gpstop's diff took out (base line 5035), and completes lane sampfront's
# CLOSED line for 4717, which its fix round outdated. Rerunnable: each edit
# says whether it was applied or already in place.
#
#   Rscript dev/rel082-backlog.R
p <- "C:/Users/adf44/source/r/frmtmb-wt-rel082/dev/test-backlog.md"
x <- sub("\r$", "", readLines(p, warn = FALSE))
n0 <- sum(grepl("CLOSED in the round of 2026-11-11", x, fixed = TRUE))
edit <- function(what, i, new) {
  if (length(i) == 1L) { x[i] <<- new; cat("applied:", what, "\n") }
  else cat("already in place (", length(i), " matches): ", what, "\n", sep = "")
}
edit("section heading",
     which(x == "## Filed in the round of 2026-11-11"),
     "## Filed at the 0.82.0 release (2026-11-11)")
stopifnot(sum(x == "## Filed at the 0.82.0 release (2026-11-11)") == 1L,
          which(x == "## Filed at the 0.82.0 release (2026-11-11)") <
            which(x == "## Filed at the 0.81.0 release (2026-11-10)"))
i <- which(x == "### MINOR and performance")
stopifnot(length(i) == 1L)
if (nzchar(x[i - 1L])) {
  x <- append(x, "", after = i - 1L); cat("applied: blank line\n")
} else cat("already in place: blank line\n")
edit("sampfront 4717",
     which(x == "CLOSED in the round of 2026-11-11 (lane sampfront): subset taken by row names."),
     paste0("CLOSED in the round of 2026-11-11 (lane sampfront): a subset is ",
            "placed by row name, or else by\n  value on the fitted row equal in ",
            "every variable the model reads."))
# the consolidation's own items, as the last subsection of the release's
# section
cons <- c(
  "### From the consolidation",
  "",
  "- **`fitted()` of a zero-inflation or hurdle predictor at its end gives NaN",
  "  standard errors with a warning** (MEDIUM, regression against 0.81.0;",
  "  merged-build conflict of lanes boundfam and ppdraw). With `zi ~ x` on",
  "  `zero_inflated_binomial()` (test-predict-draw-edge.R's `pde_zib_data(0)`)",
  "  and `hu ~ x` on `hurdle_poisson()` with no zero response, the intercept",
  "  runs to its limit and the slope is \"undefined there\". `fitted()` reads",
  "  NaN Est.Error on 195 of 200 and 297 of 300 rows and warns that they",
  "  move along a direction the fit does not determine; 0.81.0 reads 0 NaN",
  "  (median Est.Error 0.1385 and 0.2375). `se_null_drop_held()` keeps the",
  "  SEs only where every lost direction is held, and the slope is not",
  "  held. The same on lane boundfam's own build; lane ppdraw's build alone",
  "  gives 0.81.0's values. `dev/rel082-fitted-held.R`; the two warnings",
  "  escape at test-predict-draw-edge.R:101 in every tier.",
  "- **Two student() rows of lane gpstop's sweep move on lane boundfam's",
  "  build** (MEDIUM, to measure). `exp_a_student_n40_s2` falls from",
  "  10.232192 to 9.068636 (code 1 with two warnings before and after) on",
  "  the merged build and on boundfam's alone; `pow_a_student_n120_s2`",
  "  falls from 55.665287 to -4.946904 on boundfam's build alone, and the",
  "  merged build's second reading of the body brings it back to 55.633557",
  "  (0.0317 short). Every other lane's build gives 0.81.0's values.",
  "  `dev/rel082-log/kt/gpstop-out/attrib-sum.txt` (`dev/rel082-keytables.sh",
  "  gpstop-attrib`).",
  "- **frmtmb.eam and frmtmb.ode have tests that need frmtmb 0.82.0 and",
  "  floors that do not say so** (LOW). eam's test-defects.R (3 failures on",
  "  0.81.0, the limit rule's text) and ode's test-ode.R (1 error, the",
  "  `ranef()` layout) changed this round with no change to either",
  "  package's code, so neither was bumped; a check of either against",
  "  frmtmb 0.81.0 fails its tests.",
  "")
if (!any(grepl("of a zero-inflation or hurdle predictor at its end gives NaN", x, fixed = TRUE))) {
  j <- which(x == "## Filed at the 0.81.0 release (2026-11-10)")
  stopifnot(length(j) == 1L)
  x <- append(x, cons, after = j - 1L)
  cat("applied: consolidation items\n")
} else cat("already in place: consolidation items\n")
# Severities for the lanes' items that carry none (the release review,
# minor 7), by the backlog's own rule of thumb: MEDIUM where a fit or a
# reader can give a wrong or short value without a word, LOW for a
# visible refusal or warning, a parity gap, a cost or a text. Keyed by a
# piece of each item's first line; written as "**SEV. " after "- ".
sev <- c(
  "A joint runaway is still called flat" = "LOW",
  "A zero-inflation probability sent to 1 is called the family's" = "MEDIUM",
  "A negbinomial shape at its limit stays \"nocurv\" behind a" = "LOW",
  "3022 and 3037 are family fit-end checks" = "LOW",
  "A skew normal under `quadrature = TRUE` started at the stall can" = "MEDIUM",
  "The Gauss-Kronrod rule's error depends on its anchor" = "MEDIUM",
  "The nested quadrature route still gives no boundary verdict" = "LOW",
  "`sd_up_escape()` does not run under quadrature" = "MEDIUM",
  "Cross-lane (boundfam): a coefficient that runs off under one-block" = "MEDIUM",
  "The Laplace allowance still costs 1/2 per censored row on the" = "LOW",
  "hurdle_lognormal and zero_inflated_beta are outside the spike table" = "MEDIUM",
  "A free treatment contrast of `sigma ~ 1 + f` under a bounding" = "LOW",
  "An invalid prior meets the weights refusal before its own error" = "LOW",
  "The observation-level warning names a multivariate term with its" = "LOW",
  "asym_laplace with 72 of 80 rows right censored and `(1 | g4)` ends" = "MEDIUM",
  "S2. A box on a distributional parameter's own class can land one ulp" = "LOW",
  "F11 seed 8's higher maximum is at a boundary that only a refit" = "LOW",
  "The quadratic model under-reads a stop on a curved ridge tilted by a" = "LOW",
  "huber() with a group term stops short in silence on half the fits" = "MEDIUM",
  "Silent short nonlinear fits left on the review's fresh rows" = "MEDIUM",
  "`grad_headroom()` on ordinal nl fits (4793), a `frmtmb_control()`" = "LOW",
  "The two exact `gp()` rows 0.80.0's Limits name beside `student()` stay" = "MEDIUM",
  "`student()` beside an exact `gp()` has Laplace spikes" = "MEDIUM",
  "The degenerate check's \"collapsed\" test is exact-double equality," = "LOW",
  "A flat ordinal-mixture threshold with no warning at all" = "MEDIUM",
  "Ordinal mixtures with random effects run no jittered starts by" = "LOW",
  "`sar()` beside `ar(cov = FALSE)` is still refused" = "LOW",
  "`confint_varcorr()` still gives `sar()`'s rho a raw-scale" = "LOW",
  "On gaussian data a `student()` fit runs nu to about 1e9" = "LOW",
  "`skew_normal()` with `se(s, sigma = TRUE)` beside `ar(t, g)`" = "LOW",
  "`student()` with `fcor()` or `sar()` and `mi()`, plus a random effect" = "MEDIUM",
  "brms's derived `err[k]` columns are not in the draws" = "LOW",
  "`ranef()` and `coef()` of a fit keep lme4's `levels x coefs`" = "LOW",
  "`ranef()` and `coef()` of a fit answer for an exact or approximate" = "LOW",
  "`variables()` of a fit puts `b_Intercept` first under `0 +" = "LOW",
  "An exact `gp()`'s `zgp_<term>[i]` is not elementwise brms's" = "LOW",
  "`insight::get_predicted()` on `multinomial()`" = "LOW",
  "brms names an ordinal fit's `posterior_epred()` draw margin" = "LOW",
  "The `loo_*()` functions refuse a subset of the fitted rows" = "LOW",
  "`frm(drop_unused_levels = FALSE)` keeps an unused level in the model" = "LOW",
  "`log_lik(newdata =)` of draws on a subset computes every fitted" = "LOW",
  "A prediction rebuilds a group-level term's design under the" = "MEDIUM",
  "Laplace draws give `conditional_effects(method =" = "MEDIUM")
s0 <- grep("^## Filed at the 0.82.0 release", x)
s1 <- grep("^## Filed at the 0.81.0 release", x)
stopifnot(length(s0) == 1L, length(s1) == 1L)
for (k in names(sev)) {
  i <- s0 + which(startsWith(x[(s0 + 1L):(s1 - 1L)], "- ") &
                    grepl(k, x[(s0 + 1L):(s1 - 1L)], fixed = TRUE))
  if (length(i) != 1L) stop("severity key '", k, "' found ", length(i), " times")
  tag <- paste0(sev[[k]], ". ")
  if (grepl(paste0("^- \\*\\*", sev[[k]], "\\. "), x[i]) ||
        grepl(paste0("^- \\*\\*", sev[[k]], "\\.\\*\\* "), x[i])) next
  x[i] <- if (startsWith(x[i], "- **")) sub("^- \\*\\*", paste0("- **", tag), x[i])
          else sub("^- ", paste0("- **", sev[[k]], ".** "), x[i])
}
cat("severities in place:", length(sev), "\n")
# After the release review: the consolidation's three items are closed
# (lane boundfam's consolidation fixes; the eam and ode patch bumps), and
# the double stationary escape under quadrature is filed with its
# numbers (minor 6). Each CLOSED line goes after the last line of its
# item, once.
close_item <- function(key, line) {
  i <- which(startsWith(x, "- **") & grepl(key, x, fixed = TRUE))
  stopifnot(length(i) == 1L)
  if (any(x == line)) return(cat("already in place:", key, "\n"))
  j <- i
  while (j < length(x) && startsWith(x[j + 1L], "  ")) j <- j + 1L
  x <<- append(x, line, after = j)
  cat("applied: closed", key, "\n")
}
close_item("of a zero-inflation or hurdle predictor at its end gives NaN",
           paste0("CLOSED at the 0.82.0 release (lane boundfam, consolidation ",
                  "fixes): a limit direction is NaN only where the quantity ",
                  "moves along it by more than 1e-6 of its standard error; ",
                  "`zi ~ x` and `hu ~ x` read 0 NaN, no warning."))
close_item("Two student() rows of lane gpstop's sweep move on lane boundfam's",
           paste0("CLOSED at the 0.82.0 release (lane boundfam, consolidation ",
                  "fixes): the quiet limit message has the subclass ",
                  "`frmtmb_limit_fit`, which the nonlinear-start restart ",
                  "counts as a flat direction; both rows read 0.81.0's values."))
close_item("frmtmb.eam and frmtmb.ode have tests that need frmtmb 0.82.0 and",
           paste0("CLOSED at the 0.82.0 release: frmtmb.eam 0.13.2 and ",
                  "frmtmb.ode 0.7.4 need frmtmb (>= 0.82.0)."))
dbl <- c(
  "- **LOW. Under `quadrature = TRUE` the stationary escape runs twice.** It",
  "  runs on the Laplace objective before the rule is anchored (`quad_fit()`;",
  "  lane quadseam measured that order as needed, item 828) and again on the",
  "  rule (`fit_assembled()`). The second escape's record replaced the",
  "  first's; since the release review the fit keeps both",
  "  (`stationary_escape_keep()`). The cost of the second, measured by",
  "  `dev/rel082-quad-escape.R` (`dev/rel082-log/quad-escape.txt`,",
  "  `quad-escape-stall.txt`): skew_normal() with `(1 | g)`, 200 rows,",
  "  `quadrature = TRUE`; on 40 symmetric data sets the Laplace escape ran on",
  "  4 and the rule's on 0 (0 of 26.4 s of fitting); on 10 fits started at",
  "  the stall (test-quadrature-defects.R's case) the Laplace escape ran on 2",
  "  (gain 24.3 and 23.4) and the rule's on 0. The double cost was not seen;",
  "  a fit whose Laplace escape returns to alpha = 0 would pay 2 restarts on",
  "  the rule.",
  "")
if (!any(grepl("Under `quadrature = TRUE` the stationary escape runs twice", x,
               fixed = TRUE))) {
  j <- which(x == "## Filed at the 0.81.0 release (2026-11-10)")
  x <- append(x, dbl, after = j - 1L)
  cat("applied: double escape item\n")
} else cat("already in place: double escape item\n")
# The release review's short check: two LOW items, filed, behavior unchanged
short <- c(
  "- **LOW. The warning on an asym_laplace `y ~ x` fit's NaN fitted values",
  "  misdescribes the stop.** On gaussian data the fit stops at false",
  "  convergence (code 8) with every parameter \"nocurv\"; since lane",
  "  boundfam's `jc_nan_lost()` its 100 of 100 NaN fitted values warn \"...",
  "  move along a direction the fit does not determine\" (0.81.0: NaN in",
  "  silence). A warning is right, but the curvature is unusable at a stop",
  "  short of the maximum, which is not a direction the fit leaves",
  "  undetermined. That a plain asym_laplace fit stops at code 8 is older than",
  "  this round. `dev/reviews/rel082-short-probe1.R`, `short-p3.log`; release",
  "  review, short check.",
  "- **LOW. `jc_nonest()`'s move test is local to the estimate.** A `zi ~ x`",
  "  fit at zi's end (seed 7, slope -0.32) reads `predict()` at new x = -10",
  "  as 0.099 with a finite standard error 0.049, but any zi slope in about",
  "  (-7, 7) fits the data as well, since zi is about 0 on every observed",
  "  row, and a slope of -3 makes the zero-inflation probability 1 at",
  "  x = -10: the prediction there is not determined. The 1e-6 threshold",
  "  reads the derivative at the estimate only. Rows at x = -30, -60 and -100",
  "  are NaN with the warning. 0.81.0 reads the same. Release review, short",
  "  check (`dev/reviews/rel082-short-probe1.R`, `short-p3.log`).",
  "")
if (!any(grepl("misdescribes the stop.**", x, fixed = TRUE))) {
  j <- which(x == "## Filed at the 0.81.0 release (2026-11-10)")
  x <- append(x, short, after = j - 1L)
  cat("applied: short-check items\n")
} else cat("already in place: short-check items\n")
# strsplit() of "" is empty, so a blank line is kept as it is
x <- unlist(lapply(x, function(s) {
  if (grepl("\n", s, fixed = TRUE)) strsplit(s, "\n", fixed = TRUE)[[1L]] else s
}), use.names = FALSE)
stopifnot(sum(grepl("CLOSED in the round of 2026-11-11", x, fixed = TRUE)) == n0)
con <- file(p, open = "wb"); writeLines(x, con, sep = "\r\n"); close(con)
cat("CLOSED lines of the round:", n0, "\n")
