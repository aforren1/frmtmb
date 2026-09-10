# REVIEW, item 1.0b, attack 6: what the lane's 22-quantity control did
# NOT compare on an ungrouped rlddm() fit.
#
#   Rscript dev/rev-rlddm-surface.R <lib> <out.rds>
#
# Seed 4242, six learners by 60 trials, the same design
# dev/rlddm-scripts/rlddm-smoke.R uses. The lane's control compares
# numbers. This one compares what a USER SEES and what the family
# object IS, because the `ndt` link's display name changed from
# "scaled_logit(0, <ub>)" to the bare word and a printed line is a
# behaviour too.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
out <- args[[2L]]
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
    dirname(system.file(package = "frmtmb.learn")), "\n")

set.seed(4242)
d <- frm_task_design("bandit2arm", n_subject = 6, n_trial = 60,
                     seed = 4242)
s <- frm_task_simulate(rlddm(subject = id, trial = trial), d,
                       pars = list(alpha = 0.4, drift = 3, bs = 1.6,
                                   ndt = 0.2, bias = 0.5),
                       seed = 4242)[[1L]]

fit <- frm(bf(rt | dec(choice) + reward(pay1, pay2) ~ 1, drift ~ 1,
              bs ~ 1, ndt ~ 1, bias = 0.5),
           family = rlddm(subject = id, trial = trial), data = s)
fam <- frmtmb::single_response(fit)[["family"]]

txt <- function(e) tryCatch(paste(capture.output(e), collapse = "\n"),
                            error = function(err)
                              paste("ERROR:", conditionMessage(err)))
val <- function(e) tryCatch(e, error = function(err)
  paste("ERROR:", conditionMessage(err)))

rec <- list(
  logLik = as.numeric(stats::logLik(fit)),
  summary_txt = txt(print(summary(fit))),
  print_fit_txt = txt(print(fit)),
  # the family object a user can reach, before and after a fit
  link_name_fitted = val(fam[["links"]][["ndt"]][["name"]]),
  link_name_bare = val(rlddm(subject = id,
                             trial = trial)[["links"]][["ndt"]][["name"]]),
  bound_class = paste(class(fam[["ndt_bound"]]), collapse = "/"),
  bound_names = paste(names(fam[["ndt_bound"]]), collapse = ","),
  fam_slots = paste(sort(names(fam)), collapse = ","),
  # frmtmb.eam's own report on a fit from this package
  ndt_time = val(as.numeric(suppressWarnings(ndt_time(fit)))[1L]),
  # the compat table, which is user-visible documentation
  n_compat = nrow(val(frmtmb::frm_compat("rlddm"))),
  compat_status_ndtgroup =
    val(frmtmb::frm_compat("rlddm", "ndt_group()")$status),
  # the value trace and the per-trial factorization
  trace_sum = sum(log(frm_value_trace(fit)$dens)),
  # what the ndt link does at a few points, on the FITTED family
  linkinv_grid = val(as.numeric(vapply(c(-2, -1, 0, 1, 2),
                                       fam[["links"]][["ndt"]]$linkinv,
                                       numeric(1)))))
saveRDS(rec, out)
for (k in names(rec)) {
  v <- rec[[k]]
  if (is.character(v) && nchar(v) > 120) {
    cat(sprintf("%-24s <%d chars>\n", k, nchar(v)))
  } else {
    cat(sprintf("%-24s %s\n", k,
                paste(format(v, digits = 12), collapse = " ")))
  }
}
