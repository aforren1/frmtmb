# Check the INSTRUMENT. dev/sgrev-probe.R predicts which function
# UseMethod() would dispatch to by reproducing R_LookupMethod. This
# script does not predict: it CALLS the generic on an object of the
# class and reads the function that actually got the frame, out of the
# call stack at the point of the error, then compares.
#
# Usage: Rscript sgrev-dispatchcheck.R <build> <mode>
args <- commandArgs(trailingOnly = TRUE)
build <- args[[1]]; mode <- args[[2]]
LIB <- if (identical(build, "FIX")) "C:/Users/adf44/source/r/sgrev-lib" else
  "C:/Users/adf44/source/r/rellib-r3"
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
att <- function(p) suppressMessages(suppressWarnings(
  library(p, character.only = TRUE)))
OTHER7 <- c("bayesplot", "bridgesampling", "coda", "gratia", "loo",
            "posterior", "rstantools")
switch(mode,
  S = { att("brms"); att("frmtmb.sample") },
  T = { att("frmtmb.sample"); att("brms") },
  P = { for (p in OTHER7) att(p); att("frmtmb.sample") },
  Q = { att("frmtmb.sample"); for (p in OTHER7) att(p) },
  stop("mode"))

GEN <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-samplegen/dev/sgrev-out/names.rds")
who <- function(f) {
  if (!is.function(f)) return(NA_character_)
  e <- environment(f)
  if (is.null(e)) return("base")
  n <- environmentName(e); if (identical(n, "")) "<anon>" else n
}

# same lookup as the probe, copied so the two stay comparable
predict1 <- function(nm, cls) {
  g <- tryCatch(get(nm, envir = globalenv()), error = function(e) NULL)
  if (!is.function(g)) return(NA_character_)
  tbl <- tryCatch(get(".__S3MethodsTable__.", envir = environment(g),
                      inherits = FALSE), error = function(e) NULL)
  for (c1 in c(cls, "default")) {
    cand <- paste0(nm, ".", c1)
    f <- tryCatch(get(cand, envir = globalenv()), error = function(e) NULL)
    if (is.function(f)) return(paste0(who(f), "::", cand))
    if (is.environment(tbl)) {
      f <- tryCatch(get(cand, envir = tbl, inherits = FALSE),
                    error = function(e) NULL)
      if (is.function(f)) return(paste0(who(f), "::", cand))
    }
  }
  NA_character_
}

# the function that ACTUALLY received the dispatch
observe1 <- function(nm, obj) {
  g <- tryCatch(get(nm, envir = globalenv()), error = function(e) NULL)
  if (!is.function(g)) return(NA_character_)
  seen <- NA_character_
  grab <- function() {
    n <- sys.nframe()
    fs <- lapply(seq_len(n), function(i) sys.function(i))
    ix <- which(vapply(fs, function(f) identical(f, g), NA))
    if (!length(ix)) return(NA_character_)
    # min, not max: several methods here call the owner's generic
    # again from inside, so the generic appears twice in the stack and
    # max() reads the NESTED call. The frame that received THIS
    # dispatch is the one after the first occurrence.
    k <- min(ix) + 1L
    if (k > n) return(NA_character_)
    f <- fs[[k]]
    cl <- sys.call(k)
    paste0(who(f), "::", paste(deparse(cl[[1]]), collapse = ""))
  }
  # tryCatch OUTSIDE: an exiting handler established INSIDE would be
  # found first and unwind before the calling handler ever runs, which
  # is why the first spelling of this script reported 56 of 56
  # undetermined.
  tryCatch(
    withCallingHandlers(
      suppressWarnings(suppressMessages(g(obj))),
      error = function(e) {
        if (grepl("no applicable method for", conditionMessage(e),
                  fixed = TRUE) &&
            grepl(paste0("'", nm, "'"), conditionMessage(e),
                  fixed = TRUE)) {
          seen <<- "NONE::no applicable method"
        } else if (is.na(seen)) {
          seen <<- grab()
        }
      }),
    error = function(e) NULL)
  seen
}

obj_draws <- structure(list(), class = "frmtmb_draws")
obj_brms <- structure(list(), class = "brmsfit")

cat("== instrument check, dev/sgrev-dispatchcheck.R ==\n")
cat("build", build, "mode", mode, "\n")
cat(sprintf("%-20s %-6s %-46s %-46s %s\n", "generic", "class",
            "PREDICTED by the probe", "OBSERVED in the call stack",
            "agree"))
agree <- 0; disagree <- 0; undet <- 0
for (nm in GEN$names) {
  for (cl in c("frmtmb_draws", "brmsfit")) {
    obj <- if (identical(cl, "frmtmb_draws")) obj_draws else obj_brms
    p <- predict1(nm, cl)
    o <- observe1(nm, obj)
    # the observed frame is named by the call, which for a table
    # method is the bare method name
    pn <- sub("^.*::", "", p); on <- sub("^.*::", "", o)
    pp <- sub("::.*$", "", p); op <- sub("::.*$", "", o)
    if (is.na(o)) { undet <- undet + 1; next }
    if (identical(o, "NONE::no applicable method")) {
      if (is.na(p)) { agree <- agree + 1; next }
      disagree <- disagree + 1
      cat(sprintf("%-20s %-6s %-46s %-46s %s\n", nm,
                  substr(cl, 1, 6), p, o, "NO"))
      next
    }
    ok <- identical(pn, on) && identical(pp, op)
    if (ok) agree <- agree + 1 else {
      disagree <- disagree + 1
      cat(sprintf("%-20s %-6s %-46s %-46s %s\n", nm,
                  substr(cl, 1, 6), p, o, "NO"))
    }
  }
}
cat("\nagree", agree, " disagree", disagree,
    " undetermined (no error raised inside a method)", undet, "\n")
cat("DONE\n")
