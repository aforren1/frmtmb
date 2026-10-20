# Reviewer scratch (fix round 1): is running an extension's .onLoad()
# again, as local_frmtmb_registries() now does, safe for frmtmb.sample
# (frm_install_generics(), a prior-defaults provider) and for all seven?
#   Rscript dev/testiso-rev2-onload.R <scoped|plain|twice> <pkg ...>
.libPaths(c("C:/Users/adf44/source/r/wt-testiso-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
a <- commandArgs(TRUE)
mode <- a[1]
pk <- a[-1]
suppressMessages(library(frmtmb))
sys.source("C:/Users/adf44/source/r/frmtmb-wt-testiso/tests/testthat/helper-registry.R",
           envir = globalenv())
ns <- asNamespace("frmtmb")
active_before <- NULL
if (mode == "scoped") {
  f <- function() {
    local_frmtmb_registries()
    for (p in pk) loadNamespace(p)
  }
  f()
} else {
  for (p in pk) loadNamespace(p)
  if (mode == "twice") {
    # what the helper does for a namespace that appeared in its scope,
    # here without the restore in front, to see what a rerun alone adds
    for (p in pk) {
      hook <- get0(".onLoad", envir = asNamespace(p), inherits = FALSE)
      r <- tryCatch({
        hook(dirname(getNamespaceInfo(p, "path")), p); "ok"
      }, error = function(e) paste("ERROR:", conditionMessage(e)))
      cat("rerun", p, ":", r, "\n")
    }
  }
}
st <- list(
  features = length(unlist(ns$frmtmb_compat_contrib$features)),
  rules = length(ns$frmtmb_compat_contrib$rules),
  aterms = sort(names(ns$frmtmb_aterm_registry$reg)),
  checks = length(ns$frmtmb_frame_checks$fns),
  priors = length(ns$frmtmb_prior_defaults$providers))
cat("STATE", mode, paste(pk, collapse = ","), ":",
    gsub("\\s+", " ", paste(deparse(st), collapse = "")), "\n")
if ("frmtmb.sample" %in% pk) {
  sn <- asNamespace("frmtmb.sample")
  act <- Filter(function(nm) bindingIsActive(nm, sn), ls(sn, all.names = TRUE))
  cat("sample active bindings:", length(act), paste(head(act, 8), collapse = " "),
      "\n")
  # the active bindings still resolve to a function
  ok <- vapply(act, function(nm) is.function(get(nm, envir = sn)), TRUE)
  cat("all resolve to functions:", all(ok), "\n")
  set.seed(1)
  d <- data.frame(x = rnorm(40)); d$y <- 1 + d$x + rnorm(40)
  fit <- frm(bf(y ~ x) + gaussian(), data = d)
  gp <- tryCatch(nrow(get_prior(bf(y ~ x) + gaussian(), data = d)),
                 error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("get_prior rows:", gp, "\n")
}
