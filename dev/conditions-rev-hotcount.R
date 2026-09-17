# Reviewer, lane wt-conditions (round 2): how often a raise happens in
# ordinary work, counted rather than timed. frm_condition_class() runs
# once per frm_stop(), frm_warning() or frm_message(), so a trace counter
# on it counts raises, including ones a caller catches and discards.
# Also times the helper alone against stop() with a control arm.
#   Rscript dev/conditions-rev-hotcount.R      (seed 20260917)
.libPaths(c("C:/Users/adf44/source/r/conditions-lib",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.eam)})
options(frmtmb.notices = FALSE)
cnt <- new.env(); cnt$n <- 0L; cnt$who <- character()
suppressMessages(trace("frm_condition_class", where = asNamespace("frmtmb"),
  print = FALSE, tracer = quote({
    cnt$n <- cnt$n + 1L
    cnt$who <- c(cnt$who, paste(deparse(sys.call(-2L))[1L]))
  })))
set.seed(20260917)
d <- data.frame(y = rnorm(200), x = rnorm(200), g = gl(20, 10))
d$n <- 10L; d$k <- rbinom(200, 10, plogis(0.4 * d$x))
work <- list(
  gaussian_fit = quote(fg <<- frm(y ~ x, data = d)),
  binomial_fit = quote(fb <<- frm(k | trials(n) ~ x, data = d,
                                  family = binomial())),
  mixed_fit = quote(fm <<- frm(y ~ x + (1 | g), data = d)),
  predict_se = quote(predict(fm, se.fit = TRUE)),
  confint_wald = quote(confint(fm)),
  confint_profile = quote(confint(fm, method = "profile")),
  conditional_effects = quote(conditional_effects(fb, effects = "x")),
  simulate = quote(simulate(fm, nsim = 20)),
  hypothesis = quote(hypothesis(fg, "x = 0")),
  bootstrap = quote(frm_bootstrap(fg, nsim = 20, seed = 1)),
  anova_refit = quote(anova(fm, frm(y ~ 1 + (1 | g), data = d)))
)
for (nm in names(work)) {
  cnt$n <- 0L; cnt$who <- character()
  r <- tryCatch({suppressWarnings(suppressMessages(eval(work[[nm]]))); "ok"},
                error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("%-20s raises %4d  %s\n", nm, cnt$n, substr(r, 1, 60)))
  if (cnt$n) print(head(table(cnt$who), 5))
}
untrace("frm_condition_class", where = asNamespace("frmtmb"))

# the cost of the class lookup alone: interleaved, blocks past 1.2 s,
# minimum of 7 rounds, with a control arm identical to the base arm
core_env <- environment(frmtmb::frm)
eam_env <- environment(frmtmb.eam::wiener)
fcc <- frmtmb:::frm_condition_class
arm <- function(f, n) {
  t0 <- proc.time()[[3L]]
  for (i in seq_len(n)) f()
  proc.time()[[3L]] - t0
}
a_base <- function() c("error", "condition")
a_core <- function() fcc("error", core_env)
a_ext <- function() fcc("error", eam_env)
n <- 5000L
while (arm(a_core, n) < 1.2) n <- n * 2L
tab <- replicate(7L, c(base = arm(a_base, n), control = arm(a_base, n),
                       core = arm(a_core, n), ext = arm(a_ext, n)))
m <- apply(tab, 1L, min) / n * 1e6
cat("\nmicroseconds per call, n =", n, "\n")
print(round(m, 2))
cat(sprintf("control/base %.3f; lookup cost core %.2f us, extension %.2f us\n",
            m[["control"]] / m[["base"]], m[["core"]] - m[["base"]],
            m[["ext"]] - m[["base"]]))
