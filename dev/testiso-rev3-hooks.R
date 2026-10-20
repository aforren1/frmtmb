# Reviewer scratch (fix round 2): the onLoad hooks local_frmtmb_registries()
# now sets. One fresh process per case.
#   Rscript dev/testiso-rev3-hooks.R <error|other|nested|added|order>
.libPaths(c("C:/Users/adf44/source/r/wt-testiso-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
case <- commandArgs(TRUE)[1]
suppressMessages(library(frmtmb))
sys.source("C:/Users/adf44/source/r/frmtmb-wt-testiso/tests/testthat/helper-registry.R",
           envir = globalenv())
ns <- asNamespace("frmtmb")
hk <- function(p) length(getHook(packageEvent(p, "onLoad")))
allhk <- function() {
  ps <- c("frmtmb.coupling", "frmtmb.eam", "frmtmb.latent", "frmtmb.learn",
          "frmtmb.ode", "frmtmb.sample", "frmtmb.spline")
  paste(vapply(ps, hk, 1L), collapse = ",")
}
nrules <- function() length(ns$frmtmb_compat_contrib$rules)
cat("hooks before:", allhk(), "\n")
switch(case,
  error = {
    f <- function() {
      local_frmtmb_registries()
      loadNamespace("frmtmb.eam")
      stop("planted")
    }
    r <- tryCatch(f(), error = function(e) conditionMessage(e))
    cat("error caught:", r, "| dec registered:",
        "dec" %in% names(ns$frmtmb_aterm_registry$reg), "\n")
  },
  other = {
    fired <- 0L
    setHook(packageEvent("frmtmb.eam", "onLoad"),
            function(...) fired <<- fired + 1L)
    cat("hooks with another package's:", allhk(), "\n")
    f <- function() {
      local_frmtmb_registries()
      loadNamespace("frmtmb.eam")
    }
    f()
    cat("other hook fired:", fired, "\n")
  },
  nested = {
    f <- function() {
      local_frmtmb_registries()
      g <- function() {
        local_frmtmb_registries()
        loadNamespace("frmtmb.learn")
      }
      g()
      cat("after inner: hooks", allhk(), "| rules", nrules(), "\n")
      loadNamespace("frmtmb.sample")
    }
    f()
    cat("plain would give rules: eam 2 + learn 1 + sample 1 = 4\n")
  },
  added = {
    # a hook another party adds for a frmtmb package DURING the scope
    f <- function() {
      local_frmtmb_registries()
      setHook(packageEvent("frmtmb.latent", "onLoad"),
              function(...) cat("THIRD PARTY HOOK RAN\n"))
    }
    f()
    cat("third-party hook kept after scope:", hk("frmtmb.latent"), "\n")
    loadNamespace("frmtmb.latent")
  },
  order = {
    f <- function() {
      local_frmtmb_registries()
      loadNamespace("frmtmb.sample")
      loadNamespace("frmtmb.learn")
    }
    f()
    a <- ns$frmtmb_compat_contrib$features
    cat("feature order scoped:", paste(vapply(a, function(x) names(x)[1], ""),
                                      collapse = " "), "\n")
  },
  stop("unknown case"))
cat("hooks after:", allhk(), "| rules", nrules(), "| aterms",
    paste(sort(names(ns$frmtmb_aterm_registry$reg)), collapse = " "), "\n")
