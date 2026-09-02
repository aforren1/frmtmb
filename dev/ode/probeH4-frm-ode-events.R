# Probe H4: frm_ode(events = ) end to end.
#
# 1. multi-dose trajectory against the analytic superposition
# 2. gradient of a taped nl-style body against central differences
# 3. event_scale (bioavailability) as an estimated dose multiplier
# 4. infusion against the analytic infusion solution
# 5. replace / multiply
# 6. a population fit with random ka and ke, recovering the truth

suppressPackageStartupMessages({
  library(RTMB)
  library(RTMBode)
  pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-dose", quiet = TRUE)
})
source("dev/ode/pk-common.R")

say <- function(...) cat(..., "\n", sep = "")
hr <- function(t) say("\n===== ", t, " =====")

multi_dose_analytic <- function(t, ka, ke, V, amt, dose_times) {
  vapply(seq_along(t), function(i) {
    keep <- dose_times <= t[i]
    if (!any(keep)) return(0)
    u <- t[i] - dose_times[keep]
    sum(amt[keep] * ka[i] / (V[i] * (ka[i] - ke[i])) *
          (exp(-ke[i] * u) - exp(-ka[i] * u)))
  }, 0)
}

# --------------------------------------------------------------------
hr("1. multi-dose trajectory vs analytic superposition")

dose_times <- c(0, 12, 24, 36)
obs <- c(0.5, 2, 6, 11, 13, 18, 25, 30, 37, 42, 48)
d <- data.frame(id = factor(rep(c("a", "b"), each = length(obs))),
                time = rep(obs, 2))
doses <- data.frame(time = dose_times[-1], state = "depot", value = 100)

got <- frm_ode(pk_dyn, init = list(100, 0), times = d$time,
               parms = list(1, 0.2, 10), group = d$id,
               states = c("depot", "central"), output = "central",
               events = doses, atol = 1e-10, rtol = 1e-10)
ref <- multi_dose_analytic(d$time, rep(1, nrow(d)), rep(0.2, nrow(d)),
                           rep(10, nrow(d)), rep(100, 4), dose_times)
say("max abs err = ", format(max(abs(got - ref)), digits = 4))
say("max rel err = ",
    format(max(abs(got - ref) / pmax(abs(ref), 1e-8)), digits = 4))

# every dose expressed as an event, with an empty depot at t0
got2 <- frm_ode(pk_dyn, init = list(0, 0), times = d$time,
                parms = list(1, 0.2, 10), group = d$id,
                states = c("depot", "central"), output = "central",
                events = data.frame(time = dose_times, state = "depot",
                                    value = 100),
                atol = 1e-10, rtol = 1e-10)
say("dose at t0 as an event, max err = ",
    format(max(abs(got2 - ref)), digits = 4))

# per-group schedules
dd <- data.frame(time = c(dose_times[-1], 12), state = "depot",
                 value = 100,
                 group = c(rep("a", 3), "b"))
got3 <- frm_ode(pk_dyn, init = list(100, 0), times = d$time,
                parms = list(1, 0.2, 10), group = d$id,
                states = c("depot", "central"), output = "central",
                events = dd, atol = 1e-10, rtol = 1e-10)
ref3 <- c(multi_dose_analytic(obs, rep(1, length(obs)),
                              rep(0.2, length(obs)),
                              rep(10, length(obs)), rep(100, 4),
                              dose_times),
          multi_dose_analytic(obs, rep(1, length(obs)),
                              rep(0.2, length(obs)),
                              rep(10, length(obs)), rep(100, 2),
                              c(0, 12)))
say("per-group schedules, max err = ",
    format(max(abs(got3 - ref3)), digits = 4))

# a no-event group is solved by the single-solve path
say("group 'b' equals the single-dose closed form: ",
    isTRUE(all.equal(unname(got3[d$id == "b"]),
                     multi_dose_analytic(obs, rep(1, length(obs)),
                                         rep(0.2, length(obs)),
                                         rep(10, length(obs)),
                                         rep(100, 2), c(0, 12)),
                     tolerance = 1e-8)))

# --------------------------------------------------------------------
hr("2. the ADJOINT through frm_ode(events = )")

cfd <- function(f, x, h = 1e-5) {
  vapply(seq_along(x), function(j) {
    xp <- x; xp[j] <- xp[j] + h
    xm <- x; xm[j] <- xm[j] - h
    (f(xp) - f(xm)) / (2 * h)
  }, 0)
}
report <- function(label, g_ad, g_ref) {
  rel <- abs(g_ad - g_ref) / pmax(abs(g_ref), 1e-8)
  say(label)
  say("  AD  : ", paste(format(g_ad, digits = 8), collapse = "  "))
  say("  ref : ", paste(format(g_ref, digits = 8), collapse = "  "))
  say("  max rel = ", format(max(rel), digits = 4), "  VERDICT: ",
      if (max(rel) < 1e-5) "MATCHES" else "*** MISMATCH ***")
}

body_sum <- function(th) {
  "c" <- ADoverload("c")
  n <- nrow(d)
  sum(frm_ode(pk_dyn, init = list(100, 0), times = d$time,
              parms = list(rep(exp(th[1]), n), rep(exp(th[2]), n),
                           rep(exp(th[3]), n)),
              group = d$id, states = c("depot", "central"),
              output = "central", events = doses,
              atol = 1e-10, rtol = 1e-10))
}
th0 <- c(0, log(0.2), log(10))
tp <- MakeTape(body_sum, th0)
f_ana <- function(th)
  sum(multi_dose_analytic(d$time, rep(exp(th[1]), nrow(d)),
                          rep(exp(th[2]), nrow(d)),
                          rep(exp(th[3]), nrow(d)), rep(100, 4),
                          dose_times))
report("d/d(lka, lke, lV)", as.numeric(tp$jacfun()(th0)),
       cfd(f_ana, th0))

# --------------------------------------------------------------------
hr("3. event_scale: an estimated dose multiplier")

body_F <- function(th) {
  "c" <- ADoverload("c")
  n <- nrow(d)
  sum(frm_ode(pk_dyn, init = list(0, 0), times = d$time,
              parms = list(rep(exp(th[1]), n), rep(exp(th[2]), n),
                           rep(exp(th[3]), n)),
              group = d$id, states = c("depot", "central"),
              output = "central",
              events = data.frame(time = dose_times, state = "depot",
                                  value = 100),
              event_scale = rep(plogis(th[4]), n),
              atol = 1e-10, rtol = 1e-10))
}
thF <- c(th0, 0.4)
tpF <- MakeTape(body_F, thF)
f_anaF <- function(th)
  sum(multi_dose_analytic(d$time, rep(exp(th[1]), nrow(d)),
                          rep(exp(th[2]), nrow(d)),
                          rep(exp(th[3]), nrow(d)),
                          rep(plogis(th[4]) * 100, 4), dose_times))
say("value = ", format(tpF(thF), digits = 10),
    "  analytic = ", format(f_anaF(thF), digits = 10))
report("d/d(lka, lke, lV, logitF)", as.numeric(tpF$jacfun()(thF)),
       cfd(f_anaF, thF))

# --------------------------------------------------------------------
hr("4. infusion vs the analytic infusion solution")

# one-compartment IV infusion of AMOUNT over DURATION into the central
# compartment, elimination ke. A(t) = R/ke (1 - exp(-ke t)) while on.
dyn1 <- function(t, y, p) {
  "c" <- RTMB::ADoverload("c")
  list(c(-p[1] * y[1]))
}
amt <- 80; dur <- 4; ke <- 0.25
inf_analytic <- function(t) {
  R <- amt / dur
  ifelse(t <= dur, R / ke * (1 - exp(-ke * t)),
         R / ke * (1 - exp(-ke * dur)) * exp(-ke * (t - dur)))
}
tt <- c(0, 1, 2, 4, 4.5, 6, 10, 20)
gi <- frm_ode(dyn1, init = list(0), times = tt, parms = list(ke),
              events = data.frame(time = 0, state = 1, value = amt,
                                  duration = dur),
              atol = 1e-10, rtol = 1e-10)
say("max abs err = ", format(max(abs(gi - inf_analytic(tt))), digits = 4))
print(data.frame(t = tt, ode = as.numeric(gi), analytic = inf_analytic(tt)))

# repeated infusions
tt2 <- c(1, 3, 5, 9, 13, 15, 20)
gi2 <- frm_ode(dyn1, init = list(0), times = tt2, parms = list(ke),
               events = data.frame(time = c(0, 12), state = 1,
                                   value = amt, duration = dur),
               atol = 1e-10, rtol = 1e-10)
ref2 <- inf_analytic(tt2) + ifelse(tt2 > 12, inf_analytic(tt2 - 12), 0)
say("repeated infusions, max err = ",
    format(max(abs(gi2 - ref2)), digits = 4))

# infusion gradient
tpi <- MakeTape(function(th) {
  "c" <- ADoverload("c")
  sum(frm_ode(dyn1, init = list(0), times = tt2,
              parms = list(exp(th[1])),
              events = data.frame(time = c(0, 12), state = 1,
                                  value = amt, duration = dur),
              event_scale = plogis(th[2]),
              atol = 1e-10, rtol = 1e-10))
}, c(log(ke), 0.3))
f_ai <- function(th) {
  R <- plogis(th[2]) * amt / dur; k <- exp(th[1])
  one <- function(t) ifelse(t <= dur, R / k * (1 - exp(-k * t)),
                            R / k * (1 - exp(-k * dur)) * exp(-k * (t - dur)))
  sum(one(tt2) + ifelse(tt2 > 12, one(tt2 - 12), 0))
}
report("infusion, d/d(log ke, logitF)",
       as.numeric(tpi$jacfun()(c(log(ke), 0.3))),
       cfd(f_ai, c(log(ke), 0.3)))

# --------------------------------------------------------------------
hr("5. replace and multiply")

g_rep <- frm_ode(dyn1, init = list(100), times = c(1, 2, 3, 4),
                 parms = list(ke),
                 events = data.frame(time = 2, state = 1, value = 50,
                                     method = "replace"),
                 atol = 1e-10, rtol = 1e-10)
ref_rep <- c(100 * exp(-ke * 1), 100 * exp(-ke * 2),
             50 * exp(-ke * 1), 50 * exp(-ke * 2))
say("replace  max err = ", format(max(abs(g_rep - ref_rep)), digits = 4))
g_mul <- frm_ode(dyn1, init = list(100), times = c(1, 2, 3, 4),
                 parms = list(ke),
                 events = data.frame(time = 2, state = 1, value = 2,
                                     method = "multiply"),
                 atol = 1e-10, rtol = 1e-10)
ref_mul <- c(100 * exp(-ke * 1), 100 * exp(-ke * 2),
             2 * 100 * exp(-ke * 2) * exp(-ke * 1),
             2 * 100 * exp(-ke * 2) * exp(-ke * 2))
say("multiply max err = ", format(max(abs(g_mul - ref_mul)), digits = 4))

tpr <- MakeTape(function(th) {
  "c" <- ADoverload("c")
  sum(frm_ode(dyn1, init = list(100), times = c(1, 2, 3, 4),
              parms = list(exp(th[1])),
              events = data.frame(time = 2, state = 1, value = 50,
                                  method = "replace"),
              atol = 1e-10, rtol = 1e-10))
}, log(ke))
f_rep <- function(th) {
  k <- exp(th[1])
  sum(c(100 * exp(-k), 100 * exp(-2 * k), 50 * exp(-k), 50 * exp(-2 * k)))
}
report("replace, d/d(log ke) - the case deSolve events get WRONG",
       as.numeric(tpr$jacfun()(log(ke))), cfd(f_rep, log(ke)))

tpm <- MakeTape(function(th) {
  "c" <- ADoverload("c")
  sum(frm_ode(dyn1, init = list(100), times = c(1, 2, 3, 4),
              parms = list(exp(th[1])),
              events = data.frame(time = 2, state = 1, value = 2,
                                  method = "multiply"),
              atol = 1e-10, rtol = 1e-10))
}, log(ke))
f_mul <- function(th) {
  k <- exp(th[1])
  sum(c(100 * exp(-k), 100 * exp(-2 * k),
        200 * exp(-3 * k), 200 * exp(-4 * k)))
}
report("multiply, d/d(log ke)",
       as.numeric(tpm$jacfun()(log(ke))), cfd(f_mul, log(ke)))

# --------------------------------------------------------------------
hr("6. population fit: repeated dosing, random ka and ke")

set.seed(11)
n_id <- 20
dose_t <- c(0, 12, 24)
tt3 <- c(0.5, 1, 2, 4, 8, 11.9, 13, 16, 23.9, 26, 30, 36)
pd <- data.frame(id = factor(rep(seq_len(n_id), each = length(tt3))),
                 time = rep(tt3, n_id))
b_ka <- rnorm(n_id, 0, 0.30); b_ke <- rnorm(n_id, 0, 0.25)
ka_i <- exp(log(1.0) + b_ka)[as.integer(pd$id)]
ke_i <- exp(log(0.2) + b_ke)[as.integer(pd$id)]
mu <- multi_dose_analytic(pd$time, ka_i, ke_i, rep(10, nrow(pd)),
                          rep(100, 3), dose_t)
pd$conc <- mu + rnorm(nrow(pd), 0, 0.30)
pd_doses <- function()
  data.frame(time = c(12, 24), state = "depot", value = 100)

say("simulated ", nrow(pd), " rows, ", n_id, " subjects, 3 doses")
t0 <- Sys.time()
fit <- frm(
  bf(conc ~ frm_ode(pk_dyn, init = list(100, 0), times = time,
                    parms = list(exp(lka), exp(lke), exp(lV)),
                    group = id, states = c("depot", "central"),
                    output = "central", events = pd_doses),
     lka ~ 1 + (1 | id), lke ~ 1 + (1 | id), lV ~ 1, nl = TRUE) +
    gaussian(),
  data = pd, start = list(beta = c(0, log(0.25), log(8))), se = TRUE)
say("fit took ", format(as.numeric(difftime(Sys.time(), t0, units = "secs")),
                        digits = 3), "s")
print(fixef(fit))
say("truth: lka = ", format(log(1.0)), "  lke = ", format(log(0.2)),
    "  lV = ", format(log(10)))
print(VarCorr(fit))
say("truth sd: lka 0.30, lke 0.25   sigma 0.30")
say("converged: ", isTRUE(fit$optim$convergence == 0))

say("\ndone.")
