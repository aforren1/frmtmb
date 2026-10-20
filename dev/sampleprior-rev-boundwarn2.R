# Reviewer of lane sampleprior, fix round 1: two follow-ups of
# sampleprior-rev-boundwarn.R. (1) Does a class-level lb row bound a
# coefficient whose density sits on its own coef row (brms: yes, no
# warning)? (2) uniform(-1, 1) on class ar with unique time points.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-cases.R")
d <- D$wide
d$z <- rnorm(nrow(d))
pr <- suppressWarnings(set_prior("", class = "b", lb = 0) +
                         set_prior("gamma(2, 1)", class = "b", coef = "x"))
od <- cap_obj(y ~ x + z, d, gaussian(), prior = pr)
cat("par:", names(od$par), "\nlower:", CAP$args$lower, "\n")
d$t <- ave(seq_len(nrow(d)), d$g, FUN = seq_along)
w <- character(0)
p <- withCallingHandlers(set_prior("uniform(-1, 1)", class = "ar"),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
cat("frmtmb set_prior ar warnings:", length(w), "\n")
wb <- character(0)
v <- withCallingHandlers(brms::validate_prior(
  brms::set_prior("uniform(-1, 1)", class = "ar"), y ~ x + ar(time = t, gr = g),
  data = d, family = gaussian()),
  warning = function(x) { wb <<- c(wb, conditionMessage(x)); invokeRestart("muffleWarning") })
cat("brms ar warnings:", length(wb), "\n")
print(v[v$class == "ar", c("prior", "class", "lb", "ub")])
od <- cap_obj(y ~ x + ar(time = t, gr = g), d, gaussian(), prior = p)
cat("frm_sample ar: par", names(od$par), " lower", CAP$args$lower, " upper", CAP$args$upper, "\n")
