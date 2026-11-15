# Short check of lane gpfield's fix round, F2: the fallback. Where
# gp_field_newton() does not settle (forced: at most 1 step on a bernoulli
# field) or meets a curvature that is not positive (student() with
# outliers), the tape's solve decides: the draws must then be the tape
# route's bit for bit. binomial with trials() as a third family.
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
ns <- asNamespace("frmtmb")
set_steps <- function(k) {
  unlockBinding("gp_field_newton_steps", ns)
  assign("gp_field_newton_steps", k, envir = ns)
  lockBinding("gp_field_newton_steps", ns)
}
orig_newton <- ns$gp_field_newton
orig_plan <- ns$gp_field_newton_plan
tally <- new.env()
gf_put("gp_field_newton", function(fs, nt) {
  r <- orig_newton(fs, nt)
  tally$n <- (tally$n %||% 0L) + 1L
  if (is.null(r)) tally$null <- (tally$null %||% 0L) + 1L
  r
})
run <- function(tag, fit, newd, nd = 40) {
  tally$n <- 0L; tally$null <- 0L
  set.seed(3); p_def <- predict(fit, newdata = newd, summary = FALSE, ndraws = nd)
  n_def <- c(tally$n, tally$null)
  tally$n <- 0L; tally$null <- 0L
  set_steps(1L)
  set.seed(3); p_s1 <- predict(fit, newdata = newd, summary = FALSE, ndraws = nd)
  n_s1 <- c(tally$n, tally$null)
  set_steps(8L)
  gf_put("gp_field_newton_plan", function(fit) NULL)
  set.seed(3); p_tape <- predict(fit, newdata = newd, summary = FALSE, ndraws = nd)
  gf_put("gp_field_newton_plan", orig_plan)
  cat(sprintf(paste0("F2B %-10s default: newton %d calls, %d fell back | ",
                     "steps<=1: %d calls, %d fell back | steps<=1 == tape %s | ",
                     "default vs tape max |d| %.2e (sd of draws %.3g)\n"),
              tag, n_def[1], n_def[2], n_s1[1], n_s1[2],
              identical(p_s1, p_tape), max(abs(p_def - p_tape)),
              max(apply(p_tape, 2, sd))))
}
set.seed(9961)
x <- sort(round(runif(200, 0, 10), 3))
d <- data.frame(x = x, y = rbinom(200, 1, plogis(2 * sin(3 * x))))
f <- gf_fit(bf(y ~ gp(x)), data = d, family = bernoulli())$fit
run("bern200", f, data.frame(x = c(x[5], 11, 30)))
set.seed(9962)
x <- sort(round(runif(80, 0, 10), 3))
d <- data.frame(x = x, tr = 20, y = rbinom(80, 20, plogis(3 * sin(x))))
f <- gf_fit(bf(y | trials(tr) ~ gp(x)), data = d, family = binomial())$fit
run("binom80", f, data.frame(x = c(x[5], 11, 30), tr = 20))
set.seed(9963)
x <- sort(round(runif(80, 0, 10), 3))
y <- sin(x) + 0.1 * rnorm(80)
y[c(10, 40, 70)] <- y[c(10, 40, 70)] + c(3, -4, 5)
d <- data.frame(x = x, y = y)
q <- gf_fit(bf(y ~ gp(x)), data = d, family = student())
cat("student fit warn:", unique(sub(":.*", "", q$warn)), "\n")
run("student80", q$fit, data.frame(x = c(x[5], 11, 30)))
