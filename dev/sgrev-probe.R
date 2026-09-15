# Reviewer's INDEPENDENT collision probe for lane samplegen.
#
# Written from R's documented UseMethod lookup rather than from the
# lane's script: dispatch looks for `generic.class` in the environment
# the generic was CALLED from (here globalenv, so the search path),
# and then in the .__S3MethodsTable__. of environment(generic). The
# lane's probe used only the second leg. Both are reported so the
# first leg cannot hide a difference.
#
# Usage: Rscript sgrev-probe.R <build> <mode> <arm> <outfile>
#   build: FIX | BASE
#   mode:  a key in the mode table below
#   arm:   test | control   (control never loads frmtmb.sample)

args <- commandArgs(trailingOnly = TRUE)
build <- args[[1]]
mode <- args[[2]]
arm <- args[[3]]
outfile <- args[[4]]

LIB <- if (identical(build, "FIX")) {
  "C:/Users/adf44/source/r/sgrev-lib"
} else {
  "C:/Users/adf44/source/r/rellib-r3"
}
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

SAMPLE <- "frmtmb.sample"
want_sample <- identical(arm, "test")

att <- function(p) suppressMessages(suppressWarnings(
  library(p, character.only = TRUE)))
lns <- function(p) suppressMessages(suppressWarnings(
  loadNamespace(p)))

sample_att <- function() if (want_sample) att(SAMPLE)
sample_ns <- function() if (want_sample) lns(SAMPLE)

OTHER7 <- c("bayesplot", "bridgesampling", "coda", "gratia", "loo",
            "posterior", "rstantools")

notes <- character()
note <- function(...) notes <<- c(notes, paste0(...))

run_mode <- switch(
  mode,
  S = function() { att("brms"); sample_att() },
  T = function() { sample_att(); att("brms") },
  U = function() { sample_att(); lns("brms") },
  N = function() { sample_att() },
  P = function() { for (p in OTHER7) att(p); sample_att() },
  Q = function() { sample_att(); for (p in OTHER7) att(p) },
  G = function() { att("gratia"); sample_att() },
  D = function() {
    att("brms"); sample_att()
    if (want_sample) detach("package:frmtmb.sample", unload = FALSE)
    detach("package:brms", unload = FALSE)
    att("brms"); sample_att()
  },
  R = function() {
    att("loo"); sample_att()
    if ("package:loo" %in% search()) detach("package:loo", unload = FALSE)
    try(unloadNamespace("loo"), silent = TRUE)
    note("loo namespace loaded after unload: ", isNamespaceLoaded("loo"))
    att("brms")
  },
  # --- orders the lane did not run ---------------------------------
  # frmtmb.sample is never attached: only its namespace is loaded,
  # the way a package that merely Suggests it would reach it.
  RN = function() {
    att("brms")
    if (want_sample) {
      stopifnot(requireNamespace(SAMPLE, quietly = TRUE))
    }
    note("frmtmb.sample attached: ",
         "package:frmtmb.sample" %in% search())
  },
  RN2 = function() {
    if (want_sample) {
      stopifnot(requireNamespace(SAMPLE, quietly = TRUE))
    }
    att("brms")
  },
  # brms unloaded and reloaded UNDER an attached frmtmb.sample: the
  # memo in the active binding holds a namespace environment, and a
  # reload makes a new one.
  BR = function() {
    att("brms"); sample_att()
    detach("package:brms", unload = FALSE)
    try(unloadNamespace("brms"), silent = TRUE)
    note("brms namespace loaded after unload: ", isNamespaceLoaded("brms"))
    att("brms")
    note("brms reloaded, same namespace env: ", FALSE)
  },
  # core ATTACHED as well, which changes which copy of the 26
  # re-exported names is in front. Core is only Imported by
  # frmtmb.sample, so the plain orders never attach it.
  CA = function() { att("brms"); att("frmtmb"); sample_att() },
  CA2 = function() { att("brms"); sample_att(); att("frmtmb") },
  stop("unknown mode")
)

run_mode()

# --- the names, derived independently by dev/sgrev-names.R in its
# own process off the BASE install: exported by frmtmb.sample, defined
# in its namespace, body a dispatch. Not read from the lane table.
GEN <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-samplegen/dev/sgrev-out/names.rds")
base_defined <- GEN$names

# Which of those does another package own? A name is BORROWED when
# some other loaded package exports a function of the same name whose
# environment is that package's namespace.
CANDIDATES <- c("bayesplot", "bbmle", "bridgesampling", "brms", "coda",
                "gratia", "loo", "posterior", "rstan", "rstantools")

owners_of <- function(nm) {
  out <- character()
  for (p in CANDIDATES) {
    if (!isNamespaceLoaded(p)) next
    g <- tryCatch(getExportedValue(p, nm), error = function(e) NULL)
    if (!is.function(g)) next
    if (!identical(environment(g), asNamespace(p))) next
    b <- paste(deparse(body(g)), collapse = " ")
    if (grepl("UseMethod", b, fixed = TRUE)) out <- c(out, p)
  }
  out
}

# --- the lookup ---------------------------------------------------
# leg 1: a method visible from the calling environment (globalenv)
# leg 2: the method table of environment(generic)
who <- function(f) {
  if (!is.function(f)) return(NA_character_)
  e <- environment(f)
  if (is.null(e)) return("base")
  n <- environmentName(e)
  if (identical(n, "")) "<anon>" else n
}

lookup <- function(nm, cls) {
  g <- tryCatch(get(nm, envir = globalenv()), error = function(e) NULL)
  if (!is.function(g)) {
    return(list(gen = NA_character_, via = "no-generic",
                meth = NA_character_, from = NA_character_))
  }
  gen_pkg <- who(g)
  ge <- environment(g)
  tbl <- tryCatch(get(".__S3MethodsTable__.", envir = ge,
                      inherits = FALSE),
                  error = function(e) NULL)
  # R_LookupMethod tries the CALLER's environment chain and then the
  # generic's method table FOR EACH CLASS IN TURN, the real class
  # before "default". Doing all of leg 1 before leg 2 gets the order
  # wrong and scores an exported foreign `.default` as the winner over
  # a registered class method; loo exports loo_moment_match.default,
  # which is where that shows.
  for (c1 in c(cls, "default")) {
    cand <- paste0(nm, ".", c1)
    f <- tryCatch(get(cand, envir = globalenv(), inherits = FALSE),
                  error = function(e) NULL)
    if (is.function(f)) {
      return(list(gen = gen_pkg, via = "callerenv",
                  meth = c1, from = who(f)))
    }
    if (is.environment(tbl)) {
      f <- tryCatch(get(cand, envir = tbl, inherits = FALSE),
                    error = function(e) NULL)
      if (is.function(f)) {
        return(list(gen = gen_pkg, via = "s3table",
                    meth = c1, from = who(f)))
      }
    }
  }
  list(gen = gen_pkg, via = "none", meth = NA_character_,
       from = NA_character_)
}

# every class each owner has a method for, for each borrowed name
owner_classes <- function(nm, own) {
  out <- character()
  for (p in own) {
    if (!isNamespaceLoaded(p)) next
    tbl <- tryCatch(get(".__S3MethodsTable__.", envir = asNamespace(p),
                        inherits = FALSE), error = function(e) NULL)
    if (!is.environment(tbl)) next
    k <- ls(tbl, all.names = TRUE)
    k <- k[startsWith(k, paste0(nm, "."))]
    cl <- sub(paste0("^", nm, "[.]"), "", k)
    cl <- setdiff(cl, c("default", "frmtmb_draws", "frmtmb_fit"))
    out <- unique(c(out, cl))
  }
  out
}

res <- list()
for (nm in base_defined) {
  own <- owners_of(nm)
  cls <- owner_classes(nm, own)
  rows <- list()
  probe_cls <- unique(c("brmsfit", "frmtmb_draws", cls))
  for (cc in probe_cls) {
    r <- lookup(nm, cc)
    r$class <- cc
    rows[[cc]] <- r
  }
  # getS3method, checked on its own route
  gs3 <- list()
  for (cc in probe_cls) {
    f <- tryCatch(utils::getS3method(nm, cc, optional = TRUE),
                  error = function(e) NULL)
    gs3[[cc]] <- if (is.function(f)) who(f) else NA_character_
  }
  g <- tryCatch(get(nm, envir = globalenv()), error = function(e) NULL)
  res[[nm]] <- list(
    name = nm, owners = own,
    gen_pkg = who(g),
    gen_formals = if (is.function(g))
      paste(deparse(args(g)), collapse = " ") else NA_character_,
    gen_body = if (is.function(g))
      paste(deparse(body(g)), collapse = " ") else NA_character_,
    rows = rows, gs3 = gs3,
    active_ns = if (isNamespaceLoaded(SAMPLE))
      bindingIsActive(nm, asNamespace(SAMPLE)) else NA,
    active_att = if ("package:frmtmb.sample" %in% search())
      tryCatch(bindingIsActive(
        nm, as.environment("package:frmtmb.sample")),
        error = function(e) NA) else NA
  )
}

saveRDS(list(build = build, mode = mode, arm = arm,
             search = search(), notes = notes,
             base_defined = base_defined,
             base_bodies = GEN$bodies,
             res = res),
        outfile)
cat("OK", build, mode, arm, length(base_defined), "\n")
