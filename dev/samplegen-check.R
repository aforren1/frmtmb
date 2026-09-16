# Per-load-order collision probe for the 28 generics frmtmb.sample
# DEFINES. ONE PROCESS PER MODE, because the search path is built once
# per session and a single process sees one order only.
#
# The instrument is UseMethod()'s own lookup, not an error string: the
# `.__S3MethodsTable__.` of `environment(generic)`, for the class and
# then for `default`. A foreign `.default` answering is counted as a
# failure, because two of core's 27 failed silently that way.
#
#   Rscript dev/samplegen-check.R <LIB> <mode>
#
#   S  library(brms); library(frmtmb.sample)
#   T  library(frmtmb.sample); library(brms)
#   U  library(frmtmb.sample); brms only LOADED
#   N  library(frmtmb.sample); no owner loaded
#   P  every owner except brms attached, then frmtmb.sample
#   Q  frmtmb.sample, then every owner except brms attached
#   G  library(gratia); library(frmtmb.sample)
#   R  library(loo); library(frmtmb.sample); loo unloaded; brms loaded
#   I  library(brms), then a third package that IMPORTS frmtmb.sample
#      is loaded but frmtmb.sample is never attached; probed from
#      inside that package's namespace
#   D  detach and reattach both: brms, sample, detach sample, reattach,
#      detach brms, reattach
av <- commandArgs(trailingOnly = TRUE)
LIB <- av[1]
mode <- av[2]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

own28 <- c("as.mcmc", "bayes_factor", "bridge_sampler", "kfold",
           "log_lik", "log_posterior", "loo_moment_match",
           "loo_subsample", "mcmc_plot", "neff_ratio", "nsamples",
           "nuts_params", "parnames", "post_prob", "posterior_epred",
           "posterior_interval", "posterior_linpred",
           "posterior_predict", "posterior_samples", "pp_mixture",
           "predictive_error", "predictive_interval", "psis", "reloo",
           "restructure", "rhat", "stancode", "standata")
others <- c("posterior", "loo", "rstantools", "bayesplot", "coda",
            "bridgesampling", "gratia")
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
att <- function(p) quiet(library(p, character.only = TRUE))

probe_env <- globalenv()
if (mode == "S") { att("brms"); att("frmtmb.sample") }
if (mode == "T") { att("frmtmb.sample"); att("brms") }
if (mode == "U") { att("frmtmb.sample"); quiet(loadNamespace("brms")) }
if (mode == "N") { att("frmtmb.sample") }
if (mode == "P") { for (p in others) att(p); att("frmtmb.sample") }
if (mode == "Q") { att("frmtmb.sample"); for (p in others) att(p) }
if (mode == "G") { att("gratia"); att("frmtmb.sample") }
if (mode == "R") {
  att("loo"); att("frmtmb.sample")
  detach("package:loo"); unloadNamespace("loo")
  cat("LOO_UNLOADED", !isNamespaceLoaded("loo"), "\n")
  quiet(loadNamespace("brms"))
}
if (mode == "D") {
  att("brms"); att("frmtmb.sample")
  detach("package:frmtmb.sample"); att("frmtmb.sample")
  detach("package:brms"); att("brms")
}
if (mode == "I") {
  lib <- file.path(tempdir(), "samplegen-imp")
  d <- file.path(lib, "samplegenimp")
  dir.create(file.path(d, "R"), recursive = TRUE, showWarnings = FALSE)
  writeLines(c("Package: samplegenimp", "Version: 0.0.1",
               "Title: Imports frmtmb.sample", "Author: t",
               "Maintainer: t <t@t.tt>", "Description: t.",
               "License: GPL-2", "Imports: frmtmb.sample",
               "Built: R 4.6.1; ; ; windows"),
             file.path(d, "DESCRIPTION"))
  writeLines(c(sprintf("importFrom(frmtmb.sample, %s)",
                       paste(own28, collapse = ", ")),
               "export(nothing)"), file.path(d, "NAMESPACE"))
  writeLines("nothing <- function() NULL", file.path(d, "R", "samplegenimp"))
  .libPaths(c(lib, .libPaths()))
  att("brms")
  quiet(loadNamespace("samplegenimp"))
  cat("SAMPLE_ATTACHED", "package:frmtmb.sample" %in% search(), "\n")
  probe_env <- asNamespace("samplegenimp")
}

sns <- asNamespace("frmtmb.sample")
cat("MODE", mode, "\n")
cat("SAMPLE_FROM", find.package("frmtmb.sample"), "\n")
cat("CORE_FROM  ", find.package("frmtmb"), "\n")

tbl <- function(f) {
  e <- environment(f)
  if (is.null(e)) e <- baseenv()
  get0(".__S3MethodsTable__.", envir = e, inherits = FALSE)
}
gen_of <- function(g) get0(g, envir = probe_env, mode = "function")
home <- function(f) environmentName(topenv(environment(f)))

# brmsfit: is brms's class method reachable through the generic reached?
bl <- character(); bdef <- character()
for (g in own28) {
  f <- gen_of(g)
  tb <- if (is.function(f)) tbl(f) else NULL
  hit <- !is.null(tb) && exists(paste0(g, ".brmsfit"), envir = tb,
                                inherits = FALSE)
  if (!hit) {
    bl <- c(bl, g)
    if (!is.null(tb) && exists(paste0(g, ".default"), envir = tb,
                               inherits = FALSE)) bdef <- c(bdef, g)
  }
}
cat(sprintf("BRMSFIT_LOST %d of %d\n", length(bl), length(own28)))
cat("BRMSFIT_LOST_LIST", paste(bl, collapse = ","), "\n")
cat("BRMSFIT_SILENT_DEFAULT", paste(bdef, collapse = ","), "\n")

# frmtmb_draws: does the method that would run belong to frmtmb.sample?
dl <- character(); dmsg <- character()
for (g in own28) {
  f <- gen_of(g)
  mine <- get0(paste0(g, ".frmtmb_draws"), envir = sns, inherits = FALSE)
  tb <- if (is.function(f)) tbl(f) else NULL
  m <- if (is.null(tb)) NULL else
    get0(paste0(g, ".frmtmb_draws"), envir = tb, inherits = FALSE)
  if (is.null(m) && !is.null(tb)) {
    m <- get0(paste0(g, ".default"), envir = tb, inherits = FALSE)
    if (!is.null(m)) dmsg <- c(dmsg, paste0(g, "<-", home(m), ".default"))
  }
  if (is.null(m) || !identical(m, mine)) dl <- c(dl, g)
}
cat(sprintf("DRAWS_NOT_OURS %d of %d\n", length(dl), length(own28)))
cat("DRAWS_NOT_OURS_LIST", paste(dl, collapse = ","), "\n")
cat("DRAWS_FOREIGN_DEFAULT", paste(dmsg, collapse = ","), "\n")

# every OTHER class method an owner registered on its own generic,
# for every owner that is loaded: brms is one owner among eight, and a
# coda user's as.mcmc.list or loo's psis.default is lost the same way
owners_of <- list(posterior = "rhat", gratia = "posterior_samples",
  coda = "as.mcmc", bridgesampling = c("bayes_factor", "bridge_sampler",
  "post_prob"), loo = c("kfold", "loo_moment_match", "loo_subsample",
  "psis"), rstantools = c("log_lik", "nsamples", "posterior_epred",
  "posterior_interval", "posterior_linpred", "posterior_predict",
  "predictive_error", "predictive_interval"), bayesplot = c("rhat",
  "log_posterior", "neff_ratio", "nuts_params"), brms = c("mcmc_plot",
  "parnames", "posterior_samples", "pp_mixture", "reloo", "restructure",
  "stancode", "standata"))
# The CONTROL is the generic the user would reach with frmtmb.sample
# off the search path, computed in this process so both arms share one
# session. A method is DAMAGE only if the control reaches it and the
# test does not: bayesplot's and posterior's rival rhat generics already
# hide each other's methods without frmtmb.sample, and that is not
# this package's to count.
ctrl_gen <- function(g) {
  for (e in setdiff(search(), "package:frmtmb.sample")) {
    f <- get0(g, envir = as.environment(e), inherits = FALSE)
    if (is.function(f)) return(f)
  }
  NULL
}
om_n <- 0L; om_lost <- character()
for (p in names(owners_of)) {
  if (!isNamespaceLoaded(p)) next
  otb <- get0(".__S3MethodsTable__.", envir = asNamespace(p),
              inherits = FALSE)
  for (g in owners_of[[p]]) {
    og <- tryCatch(getExportedValue(p, g), error = function(e) NULL)
    if (!is.function(og) || !identical(home(og), p)) next
    ms <- grep(paste0("^", gsub(".", "[.]", g, fixed = TRUE), "[.]"),
               ls(otb, all.names = TRUE), value = TRUE)
    ms <- setdiff(ms, paste0(g, ".frmtmb_draws"))
    # as.mcmc.list is a generic of its own, not a method of as.mcmc
    ms <- grep("^as[.]mcmc[.]list[.]", ms, value = TRUE, invert = TRUE)
    f <- gen_of(g)
    rt <- if (is.function(f)) tbl(f) else NULL
    cf <- if (identical(probe_env, globalenv())) ctrl_gen(g) else NULL
    ct <- if (is.function(cf)) tbl(cf) else NULL
    for (m in ms) {
      if (is.null(ct) || !exists(m, envir = ct, inherits = FALSE)) next
      om_n <- om_n + 1L
      ok <- !is.null(rt) && exists(m, envir = rt, inherits = FALSE)
      if (!ok) om_lost <- c(om_lost, paste0(p, "::", m))
    }
  }
}
cat(sprintf("OWNER_METHODS_LOST %d of %d\n", length(om_lost), om_n))
cat("OWNER_METHODS_LOST_LIST", paste(om_lost, collapse = ","), "\n")

if (mode == "G") {
  f <- gen_of("posterior_samples")
  tb <- tbl(f)
  cat("GAM_REACHABLE", exists("posterior_samples.gam", envir = tb,
                              inherits = FALSE), "\n")
  cat("GAM_GENERIC_FROM", home(f), "\n")
}
if (mode == "R") {
  g <- gen_of("psis")
  cat("PSIS_FROM", home(g), "\n")
  cat("PSIS_STALE", !identical(g, getExportedValue("loo", "psis")), "\n")
}

rs <- vapply(own28, function(g) home(gen_of(g)), "")
tb <- table(rs)
cat("RESOLVES", paste(sprintf("%s=%d", names(tb), as.integer(tb)),
                      collapse = " "), "\n")
act <- function(nm, e) tryCatch(bindingIsActive(nm, e),
                                error = function(x) NA)
cat("ACTIVE_NS", act("log_lik", sns), "\n")
if ("package:frmtmb.sample" %in% search()) {
  cat("ACTIVE_ATTACHED",
      act("log_lik", as.environment("package:frmtmb.sample")), "\n")
}
if (mode == "I") {
  ie <- parent.env(asNamespace("samplegenimp"))
  cat("ACTIVE_IMPORTER", act("log_lik", ie), "\n")
}
cat("CHILDOK\n")
