# genrev: constructions where the adoption mechanism could CORRUPT
# rather than help.  One case per process.
args <- commandArgs(trailingOnly = TRUE)
LIB <- args[1]; case <- args[2]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(...) suppressMessages(suppressWarnings(...))

own <- function(f) {
  if (!is.function(f)) return("<absent>")
  e <- environment(f)
  if (is.null(e)) return("base")
  n <- environmentName(topenv(e)); if (nzchar(n)) n else "<anon>"
}
# does the generic `g` that a user reaches find a method for `cls`?
reach <- function(g, cls) {
  f <- tryCatch(get(g, envir = globalenv()), error = function(e) NULL)
  if (!is.function(f)) return(c("<absent>", "<absent>"))
  genv <- environment(f)
  tbl <- tryCatch(get(".__S3MethodsTable__.", envir = genv),
                  error = function(e) NULL)
  nm <- paste0(g, ".", cls)
  m <- if (!is.null(tbl) && exists(nm, envir = tbl, inherits = FALSE))
    get(nm, envir = tbl) else NULL
  c(own(f), if (is.null(m)) "<none>" else own(m))
}
# is the environment a generic's closure lives in a namespace that has
# since been unloaded and replaced?
stale <- function(g) {
  f <- tryCatch(get(g, envir = globalenv()), error = function(e) NULL)
  if (!is.function(f)) return(NA)
  e <- environment(f)
  n <- environmentName(topenv(e))
  if (!nzchar(n) || !isNamespaceLoaded(n)) return(TRUE)
  !identical(topenv(e), asNamespace(n))
}
say <- function(...) cat(sprintf(...))

if (case == "reexport") {
  # A third package re-exports frmtmb's generic (as frmtmb.sample does)
  # and an owner's namespace loads AFTERWARDS.  frmtmb's own bindings
  # migrate; does the re-exported copy?
  gens <- c("loo", "ngrps", "posterior_summary", "as_draws_df",
            "pp_check", "VarCorr", "fixef")
  q(library(genrevdep))
  say("BEFORE  frmtmb loaded=%s  brms loaded=%s\n",
      isNamespaceLoaded("frmtmb"), isNamespaceLoaded("brms"))
  for (g in gens) say("  pre  %-20s user=%-12s inside-dep=%s\n", g,
                      reach(g, "brmsfit")[1], genrevdep::genrev_owner(g))
  q(loadNamespace("brms"))
  for (g in gens) {
    r <- reach(g, "brmsfit")
    say("  post %-20s user=%-12s method=%-12s inside-dep=%s  ns-frmtmb=%s\n",
        g, r[1], r[2], genrevdep::genrev_owner(g),
        own(get(g, envir = asNamespace("frmtmb"))))
  }
}

if (case == "ownerreload") {
  # frmtmb adopts an owner at load; that owner is then unloaded and a
  # FRESH copy loaded.  No hook was left for an owner that was already
  # present, so nothing re-adopts.
  q(library(frmtmb)); q(loadNamespace("loo")); q(loadNamespace("posterior"))
  say("adopted-at-load loo=%s posterior=%s\n",
      own(get("loo", envir = asNamespace("frmtmb"))),
      own(get("ndraws", envir = asNamespace("frmtmb"))))
  ok <- tryCatch({ unloadNamespace("posterior"); unloadNamespace("loo"); TRUE },
                 error = function(e) conditionMessage(e))
  say("unloaded: %s\n", paste(ok, collapse = " "))
  q(loadNamespace("brms"))   # pulls a FRESH loo and posterior
  for (g in c("loo", "waic", "ndraws", "as_draws_df", "variables")) {
    r <- reach(g, "brmsfit")
    say("  %-16s generic-from=%-12s method=%-12s stale=%s\n",
        g, r[1], r[2], stale(g))
  }
}

if (case == "loadall") {
  # devtools/pkgload load_all on the source tree, twice, which is what
  # a developer does all day.
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    say("pkgload absent\n"); quit(save = "no")
  }
  q(library(brms))
  hookn <- function() length(getHook(packageEvent("lme4", "onLoad")))
  say("hooks-before %d\n", hookn())
  q(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-generics",
                      export_all = FALSE, quiet = TRUE))
  say("after-1 hooks=%d loo=%s postsum=%s\n", hookn(),
      reach("loo", "brmsfit")[1], reach("posterior_summary", "brmsfit")[1])
  say("after-1 method loo=%s postsum=%s\n",
      reach("loo", "brmsfit")[2], reach("posterior_summary", "brmsfit")[2])
  q(pkgload::load_all("C:/Users/adf44/source/r/frmtmb-wt-generics",
                      export_all = FALSE, quiet = TRUE))
  say("after-2 hooks=%d loo=%s postsum=%s\n", hookn(),
      reach("loo", "brmsfit")[1], reach("posterior_summary", "brmsfit")[1])
  q(pkgload::unload("frmtmb"))
  say("after-unload hooks=%d\n", hookn())
}

if (case == "defaultswap") {
  # frmtmb registers posterior_summary.default and loo_compare.default
  # in its OWN table.  After adoption those are unreachable: the
  # owner's default answers instead.  Does it answer the same?
  set.seed(1)
  m <- cbind(a = rnorm(500), b = rnorm(500, 3, 2))
  q(library(frmtmb))
  alone <- posterior_summary(m)
  alone_r <- posterior_summary(m, robust = TRUE)
  ownr <- own(get("posterior_summary", envir = globalenv()))
  q(loadNamespace("brms"))
  after <- posterior_summary(m)
  after_r <- posterior_summary(m, robust = TRUE)
  say("posterior_summary owner before=%s after=%s\n", ownr,
      own(get("posterior_summary", envir = globalenv())))
  say("cols before=%s after=%s\n", paste(colnames(alone), collapse = ","),
      paste(colnames(after), collapse = ","))
  say("identical plain=%s robust=%s\n",
      identical(alone, after), identical(alone_r, after_r))
  say("max abs diff plain=%.3e robust=%.3e\n",
      max(abs(alone - after)), max(abs(alone_r - after_r)))
  # loo_compare.default
  say("loo_compare owner=%s\n", own(get("loo_compare", envir = globalenv())))
  say("frmtmb loo_compare.default exists=%s\n",
      exists("loo_compare.default",
             envir = get(".__S3MethodsTable__.", envir = asNamespace("frmtmb")),
             inherits = FALSE))
}
say("GENREVDONE\n")
