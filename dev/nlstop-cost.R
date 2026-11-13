# Lane nlstop: the cost of the lane's two additions in one process. Arm
# "on" is the lane's build; "off" and "off2" (the control, which must read
# 1.0 against "off") replace nl_limit_message() by NULL and
# prior_conflict_escape() by the identity, which is 0.81.0's path. The
# arms rotate per round, 3 rounds; each block is the sum over its rows.
# Rows: the 450-row scan's bodies at seed 101, slope 0.6, 9 families (45
# fits), and the prior rows' linear fits at seed 6201 (21 fits).
# Usage: NLSTOP_ARM=after Rscript dev/nlstop-cost.R
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
orig <- list(lim = ns$nl_limit_message, esc = ns$prior_conflict_escape)
set_arm <- function(arm) {
  on <- arm == "on"
  utils::assignInNamespace("nl_limit_message",
    if (on) orig$lim else function(fit) NULL, "frmtmb")
  utils::assignInNamespace("prior_conflict_escape",
    if (on) orig$esc else function(obj, opt, ...) opt, "frmtmb")
}
F <- list(cumulative = cumulative(), sratio = sratio(), cratio = cratio(),
          acat = acat(), hurdle_cumulative = hurdle_cumulative(),
          bernoulli = bernoulli(), poisson = poisson(), gaussian = gaussian(),
          gamma = Gamma("log"))
bodies <- list(
  exp = c("b1 * exp(b2 * x)", "b1", "b2"),
  sat = c("b1 * (1 - exp(-b2 * x))", "b1", "b2"),
  logi = c("b1 / (1 + exp(b2 * x))", "b1", "b2"),
  prod = c("a + b1 * b2 * x", "a", "b1", "b2"),
  two = c("b1 * exp(b2 * x) + b3 * z", "b1", "b2", "b3"))
gen <- function(fam, seed, slope) {
  set.seed(seed); n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3))
  lat <- 1.5 * exp(slope * d$x) + 0.3 * d$z
  d$y <- switch(fam,
    cumulative = , sratio = , cratio = , acat =
      ordered(cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4)),
    hurdle_cumulative = ordered(ifelse(runif(n) < 0.2, 0L, as.integer(
      cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf), labels = 1:4))),
      levels = 0:4),
    bernoulli = rbinom(n, 1, plogis(lat - 2)), poisson = rpois(n, lat),
    gaussian = lat + rnorm(n, 0, 0.5), gamma = rgamma(n, 3, 3 / lat))
  d
}
jobs <- list()
for (fam in names(F)) for (bn in names(bodies)) {
  b <- bodies[[bn]]
  fm <- do.call(bf, c(list(as.formula(paste("y ~", b[1]))),
                      lapply(b[-1], function(p) as.formula(paste(p, "~ 1"))),
                      list(nl = TRUE)))
  jobs[[length(jobs) + 1L]] <- list(set = "nl", fm = fm, fam = F[[fam]],
                                    d = gen(fam, 101, 0.6), prior = NULL)
}
set.seed(6201)
dp <- data.frame(x = rnorm(200))
dp$y <- 0.5 * dp$x + rnorm(200)
for (s in c(1e2, 1e4, 1e6)) {
  dp$yl <- s * (1 + 0.05 * dp$x) + rnorm(200, 0, s * 1e-3)
  for (psd in c(1e-3, 1, 10, 30, 100, 300, 1e3)) {
    pr <- c(set_prior(sprintf("normal(0, %g)", psd), class = "b"),
            set_prior(sprintf("normal(%g, %g)", s, s), class = "Intercept"))
    jobs[[length(jobs) + 1L]] <- list(set = "prior", fm = bf(yl ~ x),
                                      fam = gaussian(), d = dp, prior = pr)
  }
}
# the prior set is repeated so that a block passes 1.2 s
run_block <- function(set) {
  t0 <- proc.time()[["elapsed"]]
  for (rep in seq_len(if (set == "prior") 8L else 1L)) for (j in jobs) {
    if (j$set != set) next
    suppressWarnings(suppressMessages(
      frm(j$fm, family = j$fam, data = j$d, prior = j$prior)))
  }
  proc.time()[["elapsed"]] - t0
}
arms <- c("off", "on", "off2")
res <- list()
for (round in 1:3) {
  ord <- arms[(seq_along(arms) + round - 2L) %% 3L + 1L]
  for (a in ord) {
    set_arm(sub("2$", "", a))
    for (set in c("nl", "prior")) {
      res[[length(res) + 1L]] <- data.frame(round = round, arm = a,
                                            set = set, t = run_block(set))
    }
  }
}
set_arm("on")
X <- do.call(rbind, res)
print(X)
for (set in c("nl", "prior")) {
  m <- tapply(X$t[X$set == set], X$arm[X$set == set], min)
  cat(sprintf("%s: min s off %.2f on %.2f off2 %.2f | on/off %.3f control off2/off %.3f\n",
              set, m[["off"]], m[["on"]], m[["off2"]], m[["on"]] / m[["off"]],
              m[["off2"]] / m[["off"]]))
}
