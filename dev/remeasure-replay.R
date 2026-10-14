# Replay one vignette of the mechanical port up to an expression and
# print the fit that expression makes: optimizer code and message, the
# warnings, the standard errors, and a finite-difference check of them.
#
#   PORT_LIB=<libs> Rscript dev/remeasure-replay.R <vignette> <id> [spell]
#
# Every expression before <id> runs as dev/brms-port/run-vignette.R runs
# it (the raw transform, or with "spell" the patches too, each seeded by
# port_seed(id)), so the data are the port's own.
HERE <- normalizePath("dev/brms-port", winslash = "/")
source(file.path(HERE, "env.R"))
source(file.path(HERE, "port-lib.R"))
a <- commandArgs(trailingOnly = TRUE)
VIG <- a[1]; TARGET <- a[2]
PATCH <- list()
if (length(a) > 2 && a[3] == "spell") source(file.path(HERE, "patches.R"))
suppressMessages(library(frmtmb))
source(file.path(HERE, "shim.R"))
cat("#", port_build(), "\n")
shim <- make_shim(suppress_brms = TRUE)
env <- new.env(parent = shim)
grDevices::pdf(NULL)
warns <- character()
msgs <- character()
val <- NULL
for (k in extract_vignette(VIG)) {
  tr <- transform_code(k$code)
  for (j in seq_along(tr)) {
    id <- sprintf("%s.%d.%d", VIG, k$idx, j)
    t <- tr[[j]]
    if (t$status %in% c("BRMS-ONLY", "PARSE-ERROR", "SETUP-SKIP")) next
    src <- if (!is.null(PATCH[[id]])) PATCH[[id]] else t$src
    set.seed(port_seed(id))
    if (identical(id, TARGET)) {
      val <- withCallingHandlers(eval(parse(text = src), envir = env),
        warning = function(w) {
          warns <<- c(warns, conditionMessage(w))
          invokeRestart("muffleWarning")
        }, message = function(m) {
          msgs <<- c(msgs, conditionMessage(m))
          invokeRestart("muffleMessage")
        })
      cat("ran:", gsub("\n *", " ", src), "\n")
      break
    }
    try(suppressWarnings(suppressMessages(
      utils::capture.output(eval(parse(text = src), envir = env)))),
      silent = TRUE)
  }
  if (!is.null(val)) break
}
stopifnot(inherits(val, "frmtmb_fit"))
cat(sprintf("logLik %.9f  convergence %d  message %s\n",
            as.numeric(logLik(val)), val$opt$convergence,
            val$opt$message %||% ""))
cat("warnings:", if (length(warns)) paste(substr(warns, 1, 200),
                                          collapse = " | ") else "none",
    "\n")
cat("messages:", if (length(msgs)) paste(substr(msgs, 1, 300),
                                         collapse = " | ") else "none",
    "\n")
se <- suppressWarnings(sqrt(diag(vcov(val))))
cat("vcov() standard errors:\n"); print(signif(se, 5))
p <- val$opt$par
cat(sprintf("max |gradient| %.3g\n", max(abs(val$obj$gr(p)))))
H <- stats::optimHess(p, val$obj$fn, val$obj$gr)
ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
cat(sprintf("finite-difference Hessian, %d parameters: eigenvalues %s\n",
            length(p), paste(signif(range(ev), 4), collapse = " to ")))
ref <- tryCatch(sqrt(diag(solve(H))), error = function(e) rep(NA, length(p)))
print(signif(stats::setNames(ref, names(p)), 5))
