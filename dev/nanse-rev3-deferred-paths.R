# Reviewer, punch round 2: does a waiting (deferred) SE report reach the
# user on ordinary paths? A natural deferral, no instrument: a gaussian
# mixed model with a 40-level factor (about 45 outer parameters, past
# se_check_np_free and se_check_share) and a (1 + x | g2) block with no
# g2 variation, which loses its SEs. For each path a fresh fit; counts
# the SE warnings that reach an outer recorder during the path, then
# during a later vcov() and summary(). "user declined" paths are meant
# to consume the report.
#   Rscript dev/nanse-rev3-deferred-paths.R lane|merge
arm <- commandArgs(TRUE)[1]
libs <- switch(arm,
  lane = c("C:/Users/adf44/source/r/wt-nanse-lib",
           "C:/Users/adf44/source/r/rellib-r5"),
  merge = c("C:/Users/adf44/source/r/nanse-rev-lib",
            "C:/Users/adf44/source/r/rellib-r6"))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("arm", arm, "from", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
phrase <- "Standard errors are not available"
mkdata <- function(s) {
  set.seed(s)
  n <- 480
  d <- data.frame(x = rnorm(n), f = factor(sample(1:40, n, TRUE)),
                  g2 = factor(rep(1:6, length.out = n)))
  d$y <- 1 + 0.5 * d$x + rnorm(40, 0, 0.5)[d$f] + rnorm(n)
  d
}
mkfit <- function(d) suppressMessages(frm(y ~ x + f + (1 + x | g2), data = d))
# pick the first seed whose fit both waits and loses an SE
for (s in 1:30) {
  d <- mkdata(s)
  w0 <- character()
  fit <- withCallingHandlers(mkfit(d), warning = function(x) {
    w0 <<- c(w0, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  if (is.null(fit$cache$se_deferred)) next
  lost <- suppressWarnings(ns$sdr_of(fit))$se_lost
  if (length(lost)) break
}
cat(sprintf("seed %d: outer pars %d, evals %s, waits, lost %s\n", s,
            length(fit$opt$par), format(fit$opt$evals),
            paste(names(lost), collapse = ",")))
count <- function(expr) {
  n <- 0L
  withCallingHandlers(expr, warning = function(x) {
    if (grepl(phrase, conditionMessage(x), fixed = TRUE)) n <<- n + 1L
    invokeRestart("muffleWarning")
  })
  n
}
user_wrapper <- function(f) suppressWarnings(fixef(f))
pkg_wrapper <- user_wrapper
environment(pkg_wrapper) <- asNamespace("stats")   # a package's function
ext_wrapper <- user_wrapper
environment(ext_wrapper) <- asNamespace("frmtmb.spline")
paths <- list(
  "fixef() at top level" = function(f) fixef(f),
  "summary() at top level" = function(f) summary(f),
  "frm_linpred(se.fit = TRUE)" = function(f) frm_linpred(f, se.fit = TRUE),
  "conditional_effects()" = function(f) conditional_effects(f, "x"),
  "user's suppressWarnings(fixef())" =
    function(f) suppressWarnings(fixef(f)),
  "global user function wrapping it" = function(f) user_wrapper(f),
  "a package's function wrapping it" = function(f) pkg_wrapper(f),
  "a frmtmb.* namespace function wrapping it" = function(f) ext_wrapper(f),
  "withCallingHandlers muffle-all" = function(f) {
    withCallingHandlers(fixef(f),
                        warning = function(w) invokeRestart("muffleWarning"))
  },
  "tryCatch(warning = ) around fixef()" = function(f) {
    r <- tryCatch(fixef(f), warning = function(w) "caught")
    if (identical(r, "caught")) message_count <<- message_count + 1L
    r
  },
  "lapply + suppressWarnings" = function(f) {
    lapply(list(f), function(z) suppressWarnings(fixef(z)))
  },
  "purrr::map" = function(f) purrr::map(list(f), fixef),
  "purrr::quietly" = function(f) {
    q <- purrr::quietly(fixef)(f)
    quiet_seen <<- sum(grepl(phrase, q$warnings, fixed = TRUE))
    q
  },
  "evaluate (knitr chunk, warning = TRUE)" = function(f) {
    out <- evaluate::evaluate("fixef(fit_e)", envir = list2env(list(fit_e = f)))
    knit_seen <<- sum(vapply(out, function(o) inherits(o, "warning") &&
                               grepl(phrase, conditionMessage(o)), TRUE))
    out
  },
  "knitr chunk warning = FALSE" = function(f) {
    e <- new.env()
    e$fit_k <- f
    knitr::knit(text = c("```{r, warning = FALSE}", "fixef(fit_k)", "```"),
                envir = e, quiet = TRUE)
  },
  "parLapply on a 1-worker cluster" = function(f) {
    cl <- parallel::makePSOCKcluster(1)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterCall(cl, function(lp) {
      .libPaths(lp); suppressMessages(library(frmtmb)); NULL
    }, .libPaths())
    parallel::parLapply(cl, list(f), function(z) frmtmb::fixef(z))
  },
  "future (multisession, 2 workers) fixef()" = function(f) {
    future::plan(future::multisession, workers = 2)
    on.exit(future::plan(future::sequential))
    v <- future::value(future::future(fixef(f), packages = "frmtmb",
                                      seed = NULL))
    v
  }
)
message_count <- 0L
quiet_seen <- NA_integer_
knit_seen <- NA_integer_
for (nm in names(paths)) {
  fit <- withCallingHandlers(mkfit(d), warning = function(x) {
    invokeRestart("muffleWarning")
  })
  stopifnot(!is.null(fit$cache$se_deferred))
  message_count <- 0L; quiet_seen <- NA_integer_; knit_seen <- NA_integer_
  a <- tryCatch(count(paths[[nm]](fit)), error = function(e) {
    cat("   error:", conditionMessage(e), "\n"); NA_integer_
  })
  pending <- !is.null(fit$cache$se_deferred)
  b <- count(vcov(fit))
  s2 <- count(summary(fit))
  extra <- c(if (message_count) "tryCatch caught it (fixef() result lost)",
             if (!is.na(quiet_seen)) paste("quietly() captured", quiet_seen),
             if (!is.na(knit_seen)) paste("chunk output shows", knit_seen))
  cat(sprintf("%-44s path %s | pending after %s | later vcov() %d, summary() %d%s\n",
              nm, format(a), pending, b, s2,
              if (length(extra)) paste0(" | ", paste(extra, collapse = "; ")) else ""))
}
