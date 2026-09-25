# Item 3.3, development check: session = against the two constructions
# it must equal. Seed 331. Output: dev/phase3b-log/session-check.txt.
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (!nzchar(Sys.getenv("P3B_INSTALLED"))) {
  suppressMessages(pkgload::load_all("extensions/frmtmb.learn",
                                     quiet = TRUE))
} else {
  library(frmtmb.learn)
}
library(frmtmb)
sink_to <- "dev/phase3b-log/session-check.txt"
out <- character(0)
say <- function(...) {
  out <<- c(out, paste0(...))
  cat(..., "\n", sep = "")
}

# 10 subjects, 2 sessions of 60 trials, trial numbers restarting
d1 <- frm_task_design("bandit2arm", n_subject = 10, n_trial = 60, seed = 331)
d2 <- frm_task_design("bandit2arm", n_subject = 10, n_trial = 60, seed = 332)
d1$session <- 1L
d2$session <- 2L
d <- rbind(d1, d2)
d$idsess <- interaction(d$id, d$session)
# draw each session from a fresh store: subject = id:session
d$choice <- frm_task_simulate(
  bandit2arm_delta(subject = idsess, trial = trial), d,
  pars = list(alpha = 0.35, tau = 3), seed = 331)[[1]]$choice

f <- bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1)
fa <- frm(f, family = bandit2arm_delta(subject = id, trial = trial,
                                       session = session), data = d)
fb <- frm(f, family = bandit2arm_delta(subject = idsess, trial = trial),
          data = d)
say("session= logLik        ", format(as.numeric(logLik(fa)), digits = 15))
say("subject=id:session     ", format(as.numeric(logLik(fb)), digits = 15))
say("identical logLik: ", identical(as.numeric(logLik(fa)),
                                    as.numeric(logLik(fb))))
say("max |coef diff|: ", format(max(abs(unlist(fixef_by_dpar(fa)) -
                                        unlist(fixef_by_dpar(fb))))))

# the objective at one parameter vector equals the two sessions' sums
p <- fa$obj$env$last.par.best
f1 <- frm(f, family = bandit2arm_delta(subject = id, trial = trial),
          data = d[d$session == 1, ])
f2 <- frm(f, family = bandit2arm_delta(subject = id, trial = trial),
          data = d[d$session == 2, ])
two <- f1$obj$fn(p) + f2$obj$fn(p)
say("objective, 2-session fit at its optimum: ",
    format(fa$obj$fn(p), digits = 15))
say("sum of the two single-session objectives at the same pars: ",
    format(two, digits = 15))
say("relative difference: ", format(abs(fa$obj$fn(p) - two) / abs(two)))

# without session =: sessions run on, so the store carries over
fc <- frm(f, family = bandit2arm_delta(subject = id, trial = trial),
          data = transform(d, trial = trial + 60 * (session - 1)))
say("no session (store carries), logLik: ",
    format(as.numeric(logLik(fc)), digits = 15))

# init at every boundary: the trace's values at each session's first
# trial are the initial values
tr <- frm_value_trace(fa)
first <- tr[tr$trial == 1, ]
say("first-trial rows: ", nrow(first), "; all q1 == 0 and q2 == 0: ",
    all(first$q1 == 0 & first$q2 == 0))
say("trace columns: ", paste(names(tr), collapse = ", "))
writeLines(out, sink_to)
