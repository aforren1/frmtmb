# One scenario per child process. REV_SCEN names it; REV_ARM picks build.
source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")

scen <- Sys.getenv("REV_SCEN")
steps <- strsplit(scen, "\\|")[[1L]]
for (s in steps) {
  if (startsWith(s, "ns:")) {
    loadNamespace(sub("^ns:", "", s))
  } else if (startsWith(s, "un:")) {
    unloadNamespace(sub("^un:", "", s))
  } else {
    suppressMessages(library(s, character.only = TRUE))
  }
}

say <- function(...) cat("P|", scen, "|", paste0(..., collapse = ""), "\n",
                        sep = "")

# --- where the visible generic comes from
vis <- tryCatch(get("log_lik"), error = function(e) NULL)
say("visible=", if (is.null(vis)) "NOTFOUND" else
      environmentName(environment(vis)))
say("rstantools_loaded=", isNamespaceLoaded("rstantools"))
if (isNamespaceLoaded("rstantools")) {
  say("is_rstantools_generic=",
      identical(vis, rstantools::log_lik))
}
say("searchpos=", paste(grep("frmtmb|brms|rstantools|loo|posterior",
                             search(), value = TRUE), collapse = ","))

# --- where are the methods registered
tbl <- function(pkg) {
  if (!isNamespaceLoaded(pkg)) return(character(0))
  e <- asNamespace(pkg)
  if (!exists(".__S3MethodsTable__.", envir = e, inherits = FALSE))
    return(character(0))
  grep("^log_lik", ls(get(".__S3MethodsTable__.", envir = e)), value = TRUE)
}
for (p in c("frmtmb", "frmtmb.sample", "rstantools", "brms")) {
  say("table.", p, "=", paste(tbl(p), collapse = ","))
}

# --- is there any .default anywhere on the path
say("any_default=",
    length(grep("^log_lik\\.default$",
                unlist(lapply(search(), function(s)
                  ls(as.environment(s)))), value = TRUE)) > 0L)

# --- a frmtmb fit must reach frmtmb's refusal
set.seed(1)
dd <- data.frame(x = rnorm(40))
dd$y <- rnorm(40, 1 + 0.5 * dd$x, 1)
fit <- frmtmb::frm(frmtmb::bf(y ~ x), family = stats::gaussian(), data = dd)
r <- tryCatch(log_lik(fit), error = function(e) e,
              warning = function(w) w)
say("fit_refusal_class=", paste(class(r), collapse = ","))
say("fit_refusal_msg=", if (inherits(r, "condition"))
      substr(conditionMessage(r), 1, 60) else
      paste("VALUE", paste(class(r), collapse = ","),
            paste(dim(r), collapse = "x")))

# --- frmtmb_draws dispatch (synthetic object: dispatch, not numerics)
if (isNamespaceLoaded("frmtmb.sample")) {
  ds <- structure(list(), class = "frmtmb_draws")
  r2 <- tryCatch(log_lik(ds), error = function(e) e)
  msg <- if (inherits(r2, "condition")) conditionMessage(r2) else "VALUE"
  say("draws_dispatch=",
      if (grepl("no applicable method", msg, fixed = TRUE)) "NO_METHOD"
      else "REACHED_METHOD")
  say("draws_msg=", substr(msg, 1, 70))
}

# --- brms's own method must still dispatch. The shipped fixtures are
# stripped of draws, so brms's own method errors on them for an unrelated
# reason; dispatch is proved by the message being IDENTICAL to the one a
# direct call to brms's method gives, and not "no applicable method".
if (isNamespaceLoaded("brms")) {
  bf1 <- get("brmsfit_example1", envir = asNamespace("brms"))
  msgof <- function(ex) tryCatch({ ex; "NOERROR" },
                                 error = function(e) conditionMessage(e))
  d1 <- msgof(brms:::log_lik.brmsfit(bf1))
  v1 <- msgof(log_lik(bf1))
  say("brms_loglik_dispatch=", identical(d1, v1),
      " noapplicable=", grepl("no applicable method", v1, fixed = TRUE))
  say("brms_loglik_msg=", substr(v1, 1, 60))
  d2 <- msgof(brms:::loo.brmsfit(bf1))
  v2 <- msgof(loo(bf1))
  say("brms_loo_dispatch=", identical(d2, v2),
      " noapplicable=", grepl("no applicable method", v2, fixed = TRUE))
}
say("DONE")
