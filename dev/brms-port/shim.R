# The data-side shim both runners evaluate the vignettes under.
#
# The eval chain is  chunk env -> shim -> globalenv -> search path.
# The shim supplies ONLY data-side conveniences (brms datasets, the one
# brms data simulator, a mirror for a dead URL). No modeling or
# post-processing function is shimmed: on the frmtmb side those must
# resolve to frmtmb or fail, which is what the audit measures.
#
# It lives in its own file so that run-vignette.R (frmtmb) and
# brms-fit.R (brms, for the estimate comparison) build the same data
# from the same code. make_shim(suppress_brms = TRUE) is the frmtmb
# side, where attaching brms would shadow every frmtmb generic.

# One seed per expression, from its id. Why: the estimate comparison
# fits frmtmb and brms in different processes, and a vignette that
# simulates its data must hand both the same rows. A seed per
# expression makes each expression's draws independent of whatever ran
# before it (sampling, a bootstrap band), so the two runs agree on
# every simulated data set and every mice() imputation.
port_seed <- function(id) {
  v <- utf8ToInt(id)
  as.integer(sum(v * seq_along(v) * 7919) %% 2147483647)
}

make_shim <- function(suppress_brms = TRUE) {
  shim <- new.env(parent = globalenv())
  for (d in c("kidney", "inhaler", "loss", "epilepsy")) {
    suppressWarnings(utils::data(list = d, package = "brms",
                                 envir = globalenv()))
  }
  shim$sim_multi_mem <- tryCatch(
    getFromNamespace("sim_multi_mem", "brms"),
    error = function(e) {
      function(nschools = 10, nstudents = 1000, change = 0.1) {
        # brms_multilevel's simulator is internal and unexported; this is
        # the documented design (two schools per student, equal weights,
        # a share of students changing school) reproduced for the audit.
        s1 <- sample(nschools, nstudents, TRUE)
        s2 <- s1
        ch <- sample(nstudents, round(change * nstudents))
        s2[ch] <- sample(nschools, length(ch), TRUE)
        eff <- stats::rnorm(nschools, 0, 3)
        y <- 20 + 0.5 * (eff[s1] + eff[s2]) + stats::rnorm(nstudents, 0, 5)
        data.frame(s1 = s1, s2 = s2, w1 = 0.5, w2 = 0.5, y = y)
      }
    }
  )
  shim$data <- function(..., package = NULL, envir = globalenv()) {
    nm <- as.character(substitute(list(...)))[-1]
    nm <- gsub('^"|"$', "", nm)
    suppressWarnings(try(utils::data(list = nm, package = package,
                                     envir = globalenv()), silent = TRUE))
    miss <- nm[!vapply(nm, exists, logical(1), envir = globalenv())]
    if (length(miss)) {
      suppressWarnings(try(utils::data(list = miss, package = "brms",
                                       envir = globalenv()), silent = TRUE))
    }
    invisible(nm)
  }
  # brms_multilevel points at a UCLA URL that no longer serves the file;
  # the same data is the author's own mirror used by brms_distreg.
  shim$read.csv <- function(file, ...) {
    if (is.character(file) && grepl("stats.idre.ucla.edu", file)) {
      file <- "https://paul-buerkner.github.io/data/fish.csv"
    }
    utils::read.csv(file, ...)
  }
  if (suppress_brms) {
    shim$library <- function(package, ...) {
      p <- tryCatch(as.character(substitute(package)), error = function(e) "")
      if (identical(p, "brms")) {
        message("[audit] library(brms) suppressed")
        return(invisible())
      }
      eval(bquote(base::library(.(as.name(p)))), envir = globalenv())
    }
    shim$require <- shim$library
  }
  shim
}
