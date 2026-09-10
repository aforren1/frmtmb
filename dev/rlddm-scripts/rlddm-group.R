# The behavioral failure, seen, and the per-group bound against it.
#
#   Rscript dev/rlddm-scripts/rlddm-group.R <lib>
#
# Seed 909. Two groups of learners whose true non-decision times are far
# apart, so that the SLOW group's truth sits above the FAST group's
# fastest response and one bound over the whole data set cannot express
# it. That is the same construction dev/ndt-scripts/ndt-before-after.R
# used for frmtmb.eam, on the family this item is about.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
  library(frmtmb.learn)
})
cat("lib  :", lib, "\n")
cat("learn:", format(packageVersion("frmtmb.learn")), "at",
    dirname(system.file(package = "frmtmb.learn")), "\n\n")

seed <- 909L
ns <- 8L
nt <- 120L
true_ndt <- rep(c(0.18, 0.40), each = ns / 2L)
d <- frm_task_design("bandit2arm", n_subject = ns, n_trial = nt,
                     seed = seed)
i <- as.integer(d$id)
s <- frm_task_simulate(rlddm(subject = id, trial = trial), d,
                       pars = list(alpha = 0.4, drift = 3, bs = 1.6,
                                   ndt = true_ndt[i], bias = 0.5),
                       seed = seed)[[1L]]
s$grp <- factor(ifelse(as.integer(s$id) <= ns / 2L, "fast", "slow"))
own <- tapply(s$rt, s$grp, min)
cat("global min(rt)      :", format(min(s$rt), digits = 9), "\n")
cat("per-group min(rt)   :", paste(names(own), format(own, digits = 9),
                                   collapse = "  "), "\n")
cat("true ndt per group  : fast 0.18  slow 0.40\n")
cat("slow truth above the GLOBAL floor:", 0.40 > min(s$rt), "\n\n")

diag1 <- function(fit) {
  dg <- frmtmb::diagnose(fit, quiet = TRUE)
  sprintf("conv=%d maxgrad=%s pdHess=%s nbadse=%d",
          dg$convergence, formatC(dg$max_grad, digits = 3, format = "g"),
          isTRUE(dg$pdHess), length(dg$bad_se))
}

# ---- arm 1: one bound over the whole data set
f_glob <- bf(rt | dec(choice) + reward(pay1, pay2) ~ 1,
             drift ~ 1, bs ~ 1, ndt ~ 0 + grp, bias = 0.5)
fit_g <- frm(f_glob, family = rlddm(subject = id, trial = trial),
             data = s, se = TRUE)
nd1 <- s[match(levels(s$grp), as.character(s$grp)), , drop = FALSE]
ndt_g <- as.numeric(suppressWarnings(
  stats::predict(fit_g, newdata = nd1, dpar = "ndt", type = "response")))
cat("--- one global bound ---\n")
cat("logLik      :", sprintf("%.6f", as.numeric(stats::logLik(fit_g))),
    "\n")
cat("fitted ndt  :", paste(levels(s$grp), sprintf("%.8f", ndt_g),
                           collapse = "  "), "\n")
cat("rel error   :", paste(levels(s$grp),
                           sprintf("%.1f%%", 100 * abs(ndt_g - c(0.18, 0.40)) /
                                     c(0.18, 0.40)), collapse = "  "), "\n")
cat("diagnose    :", diag1(fit_g), "\n\n")

# ---- arm 2: one bound per group, which is the change
f_grp <- bf(rt | dec(choice) + reward(pay1, pay2) + ndt_group(grp) ~ 1,
            drift ~ 1, bs ~ 1, ndt ~ 0 + grp, bias = 0.5)
fit_p <- frm(f_grp, family = rlddm(subject = id, trial = trial),
             data = s, se = TRUE)
ndt_p <- as.numeric(suppressWarnings(ndt_time(fit_p, newdata = nd1)))
cat("--- one bound per group ---\n")
cat("logLik      :", sprintf("%.6f", as.numeric(stats::logLik(fit_p))),
    "\n")
cat("fitted ndt  :", paste(levels(s$grp), sprintf("%.8f", ndt_p),
                           collapse = "  "), "\n")
cat("rel error   :", paste(levels(s$grp),
                           sprintf("%.1f%%", 100 * abs(ndt_p - c(0.18, 0.40)) /
                                     c(0.18, 0.40)), collapse = "  "), "\n")
cat("diagnose    :", diag1(fit_p), "\n")
cat("bound record:",
    paste(names(frmtmb::single_response(fit_p)$family$ndt_bound$floors),
          format(frmtmb::single_response(fit_p)$family$ndt_bound$floors,
                 digits = 9), collapse = "  "), "\n")
cat("sizes       :",
    paste(frmtmb::single_response(fit_p)$family$ndt_bound$sizes,
          collapse = "  "), "\n")
cat("ndt link    :", frmtmb::single_response(fit_p)$family$links$ndt$name,
    "\n")
cat("predict()   :",
    sprintf("%.8f", as.numeric(suppressWarnings(stats::predict(
      fit_p, newdata = nd1, dpar = "ndt", type = "response")))), "\n")
cat("at the wall?: fitted ndt against the group's own floor\n")
cat("  margin ms :", sprintf("%.3f", 1000 * (as.numeric(own) - ndt_p)),
    "\n\n")

# ---- the refusals and the invariances
cat("--- refusals ---\n")
say <- function(lbl, expr) {
  m <- tryCatch({ force(expr); "NO ERROR" },
                error = function(e) conditionMessage(e))
  cat(sprintf("%-26s %s\n", lbl, substr(m, 1, 120)))
}
say("max_ndt + ndt_group", frm(f_grp, family = rlddm(subject = id,
                                                     trial = trial,
                                                     max_ndt = 0.2),
                               data = s, dry_run = "frame"))
nd_unseen <- nd1
nd_unseen$grp <- factor(c("fast", "brandnew"))
say("ndt_time, unseen group", ndt_time(fit_p, newdata = nd_unseen))
say("predict, unseen group",
    stats::predict(fit_p, newdata = nd_unseen, dpar = "ndt",
                   type = "response"))
nd_nocol <- nd1
nd_nocol$grp <- NULL
say("ndt_time, column gone", ndt_time(fit_p, newdata = nd_nocol))
sb <- s
sb$pick <- sb$choice + 1L
say("ndt_group on a softmax family",
    frm(bf(pick | reward(pay1, pay2) + ndt_group(grp) ~ 1, tau ~ 1),
        family = bandit2arm_delta(subject = id, trial = trial),
        data = sb, dry_run = "frame"))
# The grouping must NOT also be a predictor for the seam's own refusal
# to be reachable: where `ndt ~ 0 + grp` puts it in a design matrix,
# frmtmb's model frame refuses a new level first and eam's "was not
# fitted to" never runs. Both are refusals; only one of them names the
# reason, and the one that names it is the one this item inherits.
f_flat <- bf(rt | dec(choice) + reward(pay1, pay2) + ndt_group(grp) ~ 1,
             drift ~ 1, bs ~ 1, ndt ~ 1, bias = 0.5)
fit_flat <- frm(f_flat, family = rlddm(subject = id, trial = trial),
                data = s)
nd_new <- nd1
nd_new$grp <- factor(c("fast", "brandnew"))
say("unseen group, ndt_time", ndt_time(fit_flat, newdata = nd_new))
say("unseen group, predict()",
    stats::predict(fit_flat, newdata = nd_new, dpar = "ndt",
                   type = "response"))
nd_gone <- nd1
nd_gone$grp <- NULL
say("column gone, ndt_time", ndt_time(fit_flat, newdata = nd_gone))
nd_na <- nd1
nd_na$grp <- factor(c("fast", NA), levels = levels(s$grp))
say("missing group, ndt_time", ndt_time(fit_flat, newdata = nd_na))
cat("in-sample ndt_time head:",
    sprintf("%.9f", utils::head(as.numeric(
      suppressWarnings(ndt_time(fit_p))), 3)), "\n")

cat("\n--- the level order does not matter ---\n")
nd_drop <- droplevels(s[s$grp == "slow", ][1:2, ])
nd_rev <- nd1
nd_rev$grp <- factor(as.character(nd_rev$grp), levels = c("slow", "fast"))
cat("subset then droplevels :",
    sprintf("%.9f", as.numeric(suppressWarnings(
      ndt_time(fit_p, newdata = nd_drop)))), "\n")
cat("levels reversed        :",
    sprintf("%.9f", as.numeric(suppressWarnings(
      ndt_time(fit_p, newdata = nd_rev)))), "\n")
s_chr <- s
s_chr$grp <- as.character(s_chr$grp)
fit_c <- frm(f_grp, family = rlddm(subject = id, trial = trial),
             data = s_chr)
cat("character spelling ll  :",
    sprintf("%.9f", as.numeric(stats::logLik(fit_c))), "  factor ll  ",
    sprintf("%.9f", as.numeric(stats::logLik(fit_p))), "\n")

cat("\n--- the bound is KEPT across a refit on fewer rows ---\n")
drop_one <- s[-which.min(s$rt), ]
fam_fit <- frmtmb::single_response(fit_p)[["family"]]
bd0 <- fam_fit[["ndt_bound"]][["floors"]]
fam2 <- fam_fit[["family_finalize"]](fam_fit, drop_one$rt,
                                     list(ndt_group = ndt_bound_key(
                                       drop_one$grp)))
cat("floors, all rows       :", format(bd0, digits = 9), "\n")
cat("floors, fastest dropped:",
    format(fam2[["ndt_bound"]][["floors"]], digits = 9), "\n")
cat("kept                   :",
    identical(bd0, fam2[["ndt_bound"]][["floors"]]), "\n")

cat("\n--- the trace and the row factorization still hold ---\n")
tr <- frm_value_trace(fit_p)
cat("sum(log(dens)) :", sprintf("%.9f", sum(log(tr$dens))), "\n")
cat("logLik(fit)    :", sprintf("%.9f", as.numeric(stats::logLik(fit_p))),
    "\n")
