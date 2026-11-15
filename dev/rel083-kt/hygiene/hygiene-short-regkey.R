# Short check of lane hygiene's fix round, F1: attacks on
# registration_key() (R/registries.R).
#   Rscript dev/reviews/hygiene-short-regkey.R <mine|base> <scratch dir>
# 0. the installed key is the lane source's
# a. cost of a closure that captured a 50 MB data frame
# b. factory inside a factory; local() inside a factory
# c. captured environments and external pointers
# d. byte-compiled against not; srcref against none
# e. load_all() twice, unload + library(), load_all() then library()
#    (extension sources copied to <scratch>, so nothing is written in
#    the lane tree)
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
scratch <- a[2]
.libPaths(c(if (arm == "mine") "C:/Users/adf44/source/r/rellib-r22" else
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, "frmtmb from", dirname(find.package("frmtmb")), "\n")
nchk <- function() length(ns$frmtmb_frame_checks$fns)
nprv <- function() length(ns$frmtmb_prior_defaults$providers)
ce <- ns$frmtmb_frame_checks
reset <- function(n) ce$fns <- ce$fns[seq_len(n)]
set.seed(1)
d <- data.frame(x = abs(stats::rnorm(40)) + 1, y = stats::rnorm(40))
fit_says <- function() tryCatch({
  suppressMessages(frm(y ~ x, data = d)); "fitted"
}, error = function(e) paste("refused:", conditionMessage(e)))

## 0. installed code is the lane's ------------------------------------------
if (arm == "mine") {
  src <- new.env()
  sys.source("C:/Users/adf44/source/r/frmtmb-wt-rel083/R/registries.R",
             envir = src, keep.source = FALSE)
  for (f in c("registration_key", "register_once")) {
    cat("INSTALLED", f, "== lane source:",
        identical(deparse(get(f, ns)), deparse(get(f, src))), "\n")
  }
}

## b. nested factories -------------------------------------------------------
n0 <- nchk()
outer <- function(mode) {
  force(mode)
  function(label) {
    force(label)
    function(spec, frame) {
      if (mode == "refuse") stop("refused by the outer value", call. = FALSE)
      invisible(NULL)
    }
  }
}
frmtmb_register_frame_check(outer("refuse")("same label"))
frmtmb_register_frame_check(outer("log")("same label"))
cat("NESTED two checks registered, added:", nchk() - n0,
    "; fit:", fit_says(), "\n")
reset(n0)
frmtmb_register_frame_check(outer("log")("same label"))
frmtmb_register_frame_check(outer("refuse")("same label"))
cat("NESTED reverse order, added:", nchk() - n0, "; fit:", fit_says(), "\n")
reset(n0)
via_local <- function(mode) {
  force(mode)
  local(function(spec, frame) {
    if (mode == "refuse") stop("refused by local()", call. = FALSE)
    invisible(NULL)
  })
}
frmtmb_register_frame_check(via_local("refuse"))
frmtmb_register_frame_check(via_local("log"))
cat("LOCAL two checks registered, added:", nchk() - n0,
    "; fit:", fit_says(), "\n")
reset(n0)

## c. captured environments and external pointers ---------------------------
mk_cfg <- function(cfg) {
  force(cfg)
  function(spec, frame) {
    if (isTRUE(cfg$refuse)) stop("refused by cfg", call. = FALSE)
    invisible(NULL)
  }
}
cfg1 <- new.env()
cfg2 <- new.env()
frmtmb_register_frame_check(mk_cfg(cfg1))
frmtmb_register_frame_check(mk_cfg(cfg2))
cfg1$refuse <- TRUE
cat("ENV two closures on two empty environments, added:", nchk() - n0,
    "; after cfg1$refuse <- TRUE, fit:", fit_says(), "\n")
reset(n0)
cfg3 <- new.env(); cfg3$refuse <- FALSE
cfg4 <- new.env(); cfg4$refuse <- TRUE
frmtmb_register_frame_check(mk_cfg(cfg3))
frmtmb_register_frame_check(mk_cfg(cfg4))
cat("ENV different contents, added:", nchk() - n0, "; fit:", fit_says(), "\n")
reset(n0)
p1 <- get("C_Cdqrls", asNamespace("stats"))$address
p2 <- get("C_DoubleCentre", asNamespace("stats"))$address
cat("XPTR classes:", class(p1), typeof(p1), "\n")
cat("XPTR two pointers identical():", identical(p1, p2),
    "; serialize equal:", identical(serialize(p1, NULL), serialize(p2, NULL)),
    "\n")
mk_ptr <- function(ptr, tag) {
  force(ptr)
  function(spec, frame) {
    if (identical(ptr, p1)) stop("refused by pointer", call. = FALSE)
    invisible(NULL)
  }
}
frmtmb_register_frame_check(mk_ptr(p1))
frmtmb_register_frame_check(mk_ptr(p2))
cat("XPTR two closures on two pointers, added:", nchk() - n0,
    "; fit:", fit_says(), "\n")
reset(n0)

## d. byte code and srcref ---------------------------------------------------
mk_plain <- function(tag) {
  force(tag)
  function(spec, frame) invisible(tag)
}
f1 <- mk_plain("a")
f2 <- compiler::cmpfun(mk_plain("a"))
cat("BYTE f2 is compiled:",
    !is.null(tryCatch(compiler::disassemble(f2), error = function(e) NULL)),
    "\n")
frmtmb_register_frame_check(f1)
frmtmb_register_frame_check(f2)
cat("BYTE same closure plain then compiled, added:", nchk() - n0, "\n")
reset(n0)
mk_helper <- function(h) {
  force(h)
  function(spec, frame) h(frame)
}
h <- function(x) invisible(NULL)
frmtmb_register_frame_check(mk_helper(h))
frmtmb_register_frame_check(mk_helper(compiler::cmpfun(h)))
cat("BYTE captured helper plain then compiled, added:", nchk() - n0, "\n")
reset(n0)
txt <- "function(spec, frame) invisible(NULL)"
s1 <- eval(parse(text = txt, keep.source = TRUE))
s2 <- eval(parse(text = txt, keep.source = FALSE))
frmtmb_register_frame_check(s1)
frmtmb_register_frame_check(s2)
cat("SRCREF same code with and without srcref, added:", nchk() - n0, "\n")
reset(n0)
frmtmb_register_frame_check(mk_helper(eval(parse(
  text = "function(x) invisible(NULL)", keep.source = TRUE))))
frmtmb_register_frame_check(mk_helper(eval(parse(
  text = "function(x) invisible(NULL)", keep.source = FALSE))))
cat("SRCREF captured helper with and without srcref, added:", nchk() - n0,
    "\n")
reset(n0)

## a. a closure that captured 50 MB ------------------------------------------
peak <- function() sum(gc()[, 6L])
timed <- function(label, expr) {
  invisible(gc(reset = TRUE)); m0 <- peak()
  t0 <- proc.time()[["elapsed"]]
  force(expr)
  t1 <- proc.time()[["elapsed"]]
  cat(sprintf("COST %-52s %.3f s, peak +%.0f MB\n", label, t1 - t0,
              peak() - m0))
}
big <- data.frame(x = stats::rnorm(6.5e6))
cat(sprintf("COST captured object %.1f MB\n",
            as.numeric(object.size(big)) / 2^20))
mk_big <- function(dat) {
  force(dat)
  function(spec, frame) invisible(NULL)
}
timed("register the 50 MB closure", frmtmb_register_frame_check(mk_big(big)))
timed("register it again (same object)",
      frmtmb_register_frame_check(mk_big(big)))
timed("register a small closure after it",
      frmtmb_register_frame_check(mk_plain("small")))
timed("10 small closures after it",
      for (i in 1:10) frmtmb_register_frame_check(mk_plain(paste0("s", i))))
timed("library(frmtmb.eam, frmtmb.ode, frmtmb.sample)",
      for (p in c("frmtmb.eam", "frmtmb.ode", "frmtmb.sample")) {
        suppressMessages(library(p, character.only = TRUE))
      })
cat("COST checks now", nchk(), "providers", nprv(), "\n")
timed("a fit with the 50 MB closure registered", fit_says())
# keep only the extensions' own checks (named namespace functions)
ce$fns <- Filter(function(f) isNamespace(environment(f)), ce$fns)
n1 <- nchk()
p1n <- nprv()
cat("COST after removing the big one: checks", n1, "\n")

## e. reloads ----------------------------------------------------------------
cp <- file.path(scratch, paste0("regkey-src-", arm))
dir.create(cp, showWarnings = FALSE)
for (p in c("frmtmb.ode", "frmtmb.sample")) {
  file.copy(file.path("C:/Users/adf44/source/r/frmtmb-wt-rel083/extensions",
                      p), cp, recursive = TRUE)
}
cat("RELOAD start: checks", nchk(), "providers", nprv(), "\n")
for (i in 1:2) {
  suppressMessages(pkgload::load_all(file.path(cp, "frmtmb.ode"),
                                     quiet = TRUE, export_all = FALSE))
  cat("RELOAD load_all(ode) #", i, ": checks", nchk(), "\n")
}
pkgload::unload("frmtmb.ode")
suppressMessages(library(frmtmb.ode))
cat("RELOAD ode unload + library(): checks", nchk(), "\n")
unloadNamespace("frmtmb.ode")
suppressMessages(library(frmtmb.ode))
cat("RELOAD ode unloadNamespace + library(): checks", nchk(), "\n")
for (i in 1:2) {
  suppressMessages(pkgload::load_all(file.path(cp, "frmtmb.sample"),
                                     quiet = TRUE, export_all = FALSE))
  cat("RELOAD load_all(sample) #", i, ": providers", nprv(), "\n")
}
pkgload::unload("frmtmb.sample")
suppressMessages(library(frmtmb.sample))
cat("RELOAD sample unload + library() after load_all: providers", nprv(),
    "\n")
unloadNamespace("frmtmb.sample")
suppressMessages(library(frmtmb.sample))
cat("RELOAD sample unloadNamespace + library(): providers", nprv(), "\n")
cat("RELOAD end: checks", nchk(), "(start", n1, ") providers", nprv(),
    "(start", p1n, ")\n")
