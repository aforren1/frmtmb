## frmtmb's names for every case of dev/brmsnames-1024-cases.R, on one
## arm, for comparison with brms's (dev/brmsnames-1024-brms.R).
##
##   Rscript dev/brmsnames-1024-frm.R base|lane <case number> [...]
##
## base reads rellib-r12 (frmtmb 0.74.0); lane reads the lane library
## first. Writes dev/brmsnames-1024-out/frm-<arm>-<n>.rds. The prior
## check runs every coefficient row of brms's default_prior() through
## frmtmb's validate_prior() (the check frm() applies) with brms's
## spelling, and hypothesis() is called on brms's variable names, split
## into class and name as brms's hypothesis() takes them. The draws are
## frmtmb.sample's frm_sample(), one chain of 40 iterations, seed 1.
arm <- commandArgs(trailingOnly = TRUE)[1L]
ids <- as.integer(commandArgs(trailingOnly = TRUE)[-1L])
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
q(library(frmtmb.sample))
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
source("dev/brmsnames-1024-cases.R")
d <- brmsnames_data()
cases <- brmsnames_cases(frmtmb::bf, frmtmb::mvbf, frmtmb::set_rescor)

try_names <- function(e) {
  tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
}
fam_of <- function(s) get(s, envir = asNamespace("frmtmb"))()

for (i in ids) {
  cs <- cases[[i]]
  cat("case", i, cs$id, "\n")
  ref <- readRDS(sprintf("dev/brmsnames-1024-out/brms-%02d.rds", i))
  fam <- fam_of(cs$family)
  out <- list(id = cs$id, frmtmb = format(packageVersion("frmtmb")))
  fit <- tryCatch(q(frm(cs$form, data = d, family = fam)),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    out$error <- conditionMessage(fit)
    saveRDS(out, sprintf("dev/brmsnames-1024-out/frm-%s-%02d.rds", arm, i))
    next
  }
  out$internal <- names(fit$frame$par_template$beta)
  out$variables <- try_names(variables(fit))
  out$fixef <- try_names(rownames(fixef(fit)))
  out$vcov <- try_names(rownames(stats::vcov(fit)))
  out$summary <- try_names(rownames(summary(fit)$fixed))
  out$confint <- try_names(rownames(q(stats::confint(fit))))
  out$coef <- try_names(lapply(stats::coef(fit), function(a) {
    if (length(dim(a)) == 3L) dimnames(a)[[3]] else colnames(a)
  }))
  out$ranef <- try_names(lapply(ranef(fit), function(a) {
    if (length(dim(a)) == 3L) dimnames(a)[[3]] else colnames(a)
  }))
  out$varcorr <- try_names(lapply(VarCorr(fit), function(a) {
    rownames(a$sd)
  }))
  out$ce <- try_names(names(q(conditional_effects(fit))))
  pr <- as.data.frame(q(default_prior(cs$form, data = d, family = fam)))
  out$prior <- pr[, c("class", "coef", "group", "resp", "dpar", "nlpar")]
  rows <- which(nzchar(ref$prior$coef))
  out$prior_ok <- vapply(rows, function(r) {
    a <- brmsnames_prior_args(ref$prior[r, ])
    p <- tryCatch(do.call(frmtmb::set_prior, a), error = function(e) e)
    if (inherits(p, "error")) return(conditionMessage(p))
    v <- tryCatch(q(validate_prior(p, cs$form, data = d, family = fam)),
                  error = function(e) e)
    if (inherits(v, "error")) return(conditionMessage(v))
    v <- as.data.frame(v)
    if (!any(v$source == "user" & v$coef == a$coef)) return("not applied")
    "ok"
  }, "")
  hv <- grep("^(b|bsp|bs|bcs|sds|sd)_", ref$variables, value = TRUE)
  out$hyp <- vapply(hv, function(v) {
    cls <- sub("_.*$", "", v)
    rest <- sub("^[a-z]+_", "", v)
    h <- tryCatch(q(hypothesis(fit, paste0(rest, " = 0"), class = cls)),
                  error = function(e) e)
    if (inherits(h, "error")) conditionMessage(h) else "ok"
  }, "")
  dr <- tryCatch(q(frm_sample(fit, chains = 1, iter = 40, warmup = 20,
                              refresh = 0, seed = 1)),
                 error = function(e) e)
  out$draws <- if (inherits(dr, "error")) {
    paste("ERROR:", conditionMessage(dr))
  } else try_names(variables(dr))
  saveRDS(out, sprintf("dev/brmsnames-1024-out/frm-%s-%02d.rds", arm, i))
}
