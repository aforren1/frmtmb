# Reviewer 2, item 8: the session-run refusal (false alarms and misses)
# and the by-name refusals in frmtmb.eam. dry_run only.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.learn); library(frmtmb.eam)})
msg <- function(expr) tryCatch({ expr; "ACCEPTED" }, error = function(e)
  paste0("REFUSED: ", substr(conditionMessage(e), 1, 140)))
f <- bf(choice | reward(pay1, pay2) ~ 1, tau ~ 1)
fam <- bandit2arm_delta(subject = id, trial = trial, session = session)
d0 <- frm_task_design("bandit2arm", n_subject = 3L, n_trial = 90L, seed = 5L)
d0$choice <- frm_task_simulate(bandit2arm_delta(subject = id, trial = trial), d0,
  pars = list(alpha = 0.4, tau = 3), seed = 5)[[1L]]$choice
run <- function(lab, sess, d = d0) {
  d$session <- sess
  cat(sprintf("%-58s %s\n", lab, msg(frm(f, family = fam, data = d, dry_run = "objective"))))
}
blk <- (d0$trial - 1L) %/% 30L + 1L
run("a b a, continuous numbering", c("a", "b", "a")[blk])
run("a a b, continuous numbering", c("a", "a", "b")[blk])
set.seed(1); sh <- sample(nrow(d0))
d_sh <- d0[sh, ]
run("a b a, rows shuffled", c("a", "b", "a")[((d_sh$trial - 1L) %/% 30L + 1L)], d_sh)
run("a a b, rows shuffled", c("a", "a", "b")[((d_sh$trial - 1L) %/% 30L + 1L)], d_sh)
run("interleaved a b a b ... (two concurrent stores)", c("a", "b")[(d0$trial %% 2L) + 1L])
# numbering restarted per session: a (1..30), b (1..30), a again (31..60)
d_r <- d0; d_r$trial <- c(1:30, 1:30, 31:60)[d0$trial]
run("a b a, restarted numbering (b restarts, a continues)", c("a", "b", "a")[blk], d_r)
# only one subject violates
s_one <- ifelse(d0$id == levels(factor(d0$id))[2], c("a", "b", "a")[blk], c("a", "a", "b")[blk])
run("violation in one subject only", s_one)

cat("\n-- frmtmb.eam refusals by name\n")
set.seed(2)
dl <- lba_simulate(200, v = c(2, 1.2), A = 0.5, k = 0.5, ndt = 0.2)
dr <- rdm_simulate(200, v = c(3, 2), A = 0.5, k = 0.5, ndt = 0.2)
dl$code <- 0L; dr$code <- 0L
chk <- function(lab, expr) cat(sprintf("%-40s %s\n", lab, msg(expr)))
chk("lba(2, contaminant = TRUE)", lba(2, contaminant = TRUE))
chk("rdm(2, contaminant = TRUE)", rdm(2, contaminant = TRUE))
chk("lba(2, contaminant = NA)", lba(2, contaminant = NA))
chk("lba + cens()", frm(rt | vint(choice) + cens(code) ~ 1, family = lba(2), data = dl,
                        dry_run = "objective"))
chk("lba + trunc(ub = 5)", frm(rt | vint(choice) + trunc(ub = 5) ~ 1, family = lba(2),
                               data = dl, dry_run = "objective"))
chk("rdm + cens() (right only)", frm(rt | vint(choice) + cens(code) ~ 1, family = rdm(2),
                                     data = dr, dry_run = "objective"))
dg <- ddm_simulate(200, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5); dg$code <- 0L
chk("gddm + cens()", frm(rt | dec(upper) + cens(code) ~ 1, family = gddm(), data = dg,
                         dry_run = "objective"))
chk("gddm + trunc(ub = 5)", frm(rt | dec(upper) + trunc(ub = 5) ~ 1, family = gddm(),
                                data = dg, dry_run = "objective"))
chk("wiener(variability = 'sv') + cens()", frm(rt | dec(upper) + cens(code) ~ 1,
    family = wiener(variability = "sv"), data = dg, dry_run = "objective"))
chk("wiener_gng + cens()", frm(rt | dec(upper) + cens(code) ~ 1, family = wiener_gng(),
                               data = dg, dry_run = "objective"))
chk("wiener(contaminant = TRUE, allow_unreachable)", wiener(contaminant = TRUE,
    allow_unreachable = TRUE))
chk("wiener(contaminant_range = c(0,1)) alone", wiener(contaminant_range = c(0, 1)))
chk("wiener(contaminant = TRUE, range c(1, 1))", wiener(contaminant = TRUE,
    contaminant_range = c(1, 1)))
