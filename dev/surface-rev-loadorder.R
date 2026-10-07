# Reviewer of lane surface, claim 1: load order and dispatch.
#   Rscript dev/surface-rev-loadorder.R <arm> <scenario>
# scenario: a "+"-separated load sequence, e.g. "frmtmb+brms",
# "brms+frmtmb.sample", "frmtmb.sample+ns:brms" (ns: = loadNamespace)
a <- commandArgs(TRUE)
source("dev/surface-rev-env.R")
rev_env(a[1])
steps <- strsplit(a[2], "+", fixed = TRUE)[[1]]
cat("arm", a[1], "| scenario", a[2], "\n")
for (s in steps) {
  msgs <- character()
  withCallingHandlers({
    if (startsWith(s, "ns:")) loadNamespace(sub("ns:", "", s))
    else library(s, character.only = TRUE)
  }, message = function(m) {
    msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
  }, packageStartupMessage = function(m) {
    msgs <<- c(msgs, paste("[startup]", conditionMessage(m)))
    invokeRestart("muffleMessage")
  })
  cat("load", s, ": ", length(msgs), "message(s)\n")
  for (x in msgs) cat("  |", gsub("\n", "\n  | ", trimws(x)), "\n")
}
cat("frmtmb", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")

where <- function(nm) {
  if (!exists(nm, envir = globalenv())) return("<not on search path>")
  f <- get(nm, envir = globalenv())
  e <- environment(f)
  paste0(environmentName(e), " (found in ",
         environmentName(as.environment(find(nm)[1])), ")")
}
for (nm in c("stancode", "standata", "pp_mixture", "add_criterion",
             "parnames", "loo", "plot")) {
  cat(sprintf("%-14s generic from %s\n", nm, where(nm)))
}

set.seed(1)
d <- data.frame(x = rnorm(40)); d$y <- d$x + rnorm(40)
f <- frmtmb::frm(frmtmb::bf(y ~ x), family = gaussian(), data = d)
set.seed(4)
dm <- data.frame(y = c(rnorm(60, -2), rnorm(60, 3)))
fm <- frmtmb::frm(frmtmb::bf(y ~ 1),
                  family = frmtmb::mixture(gaussian(), gaussian()),
                  data = dm)
call_here <- function(nm, ...) {
  if (!exists(nm, envir = globalenv())) return(invisible(NULL))
  g <- get(nm, envir = globalenv())
  rev_show(paste0(nm, "(", paste(...names(), collapse = ","), ")"),
           g(...))
}
rev_show("stancode(fit)", stancode(f))
rev_show("standata(fit)", standata(f))
p <- rev_show("pp_mixture(mixfit)", pp_mixture(fm))
if (is.array(p)) cat("   dim", dim(p), "|", dimnames(p)[[3]], "\n")
rev_show("add_criterion(fit, 'loo')", add_criterion(f, "loo"))
if ("brms" %in% loadedNamespaces()) {
  rev_show("brms::stancode(fit)", brms::stancode(f))
  rev_show("brms::pp_mixture(mixfit)", brms::pp_mixture(fm))
  rev_show("brms::add_criterion(fit,'loo')", brms::add_criterion(f, "loo"))
  # brms's own objects must still reach brms's methods
  bfit <- suppressMessages(brms::brm(y ~ x, data = d, empty = TRUE))
  sc <- rev_show("stancode(empty brmsfit)", stancode(bfit))
  rev_show("stancode(formula, data) [brms default]",
           stancode(y ~ x, data = d))
  sd <- rev_show("standata(formula, data) [brms default]",
                 standata(y ~ x, data = d))
  if (is.list(sd)) cat("   standata names:", head(names(sd), 6), "\n")
  rev_show("pp_mixture(empty brmsfit)", pp_mixture(bfit))
  rev_show("add_criterion(empty brmsfit,'loo')",
           add_criterion(bfit, "loo"))
}
if ("frmtmb.sample" %in% loadedNamespaces()) {
  ds <- suppressWarnings(frmtmb.sample::frm_sample(
    frmtmb::bf(y ~ x), family = gaussian(), data = d, chains = 1,
    iter = 300, refresh = 0, seed = 1))
  rev_show("stancode(draws)", stancode(ds))
  rev_show("standata(draws)", standata(ds))
  ds2 <- rev_show("add_criterion(draws,'loo')", add_criterion(ds, "loo"))
  if (inherits(ds2, "frmtmb_draws")) {
    cat("   stored:", names(ds2$criteria), "| loo(ds2) identical stored:",
        identical(loo(ds2), ds2$criteria$loo), "\n")
  }
  rev_show("pp_mixture(non-mixture draws)", pp_mixture(ds))
  if ("brms" %in% loadedNamespaces()) {
    rev_show("brms::stancode(draws)", brms::stancode(ds))
    rev_show("brms::add_criterion(draws,'waic')",
             brms::add_criterion(ds, "waic"))
  }
}
