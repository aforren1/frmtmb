# Review of lane hygiene, item 3154 (register_once()).
# 1. Two frame checks made by one factory (same code, different captured
#    value) are two checks. Does the lane's key collapse them?
# 2. The same for prior-default providers.
# 3. library() / unloadNamespace() / library() of frmtmb.ode.
# 4. Cost: time and memory of a registration whose closure's topenv() is
#    a namespace (every binding there is get()-forced).
#
# Usage: Rscript dev/reviews/hygiene-review-regonce.R <mine|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "mine") "C:/Users/adf44/source/r/rellib-r22" else
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
t0 <- proc.time()[["elapsed"]]
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", dirname(find.package("frmtmb")),
    "load s", round(proc.time()[["elapsed"]] - t0, 2), "\n")
ns <- asNamespace("frmtmb")
nchk <- function() length(ns$frmtmb_frame_checks$fns)
nprv <- function() length(ns$frmtmb_prior_defaults$providers)

## 1. factory checks ------------------------------------------------------
ran <- new.env()
make_check <- function(tag) {
  force(tag)
  function(spec, frame) {
    assign(tag, TRUE, envir = ran)
    # stands for any data condition an extension refuses
    if (tag == "refuse_negative_x") {
      stop("refused by check ", tag, call. = FALSE)
    }
    invisible(NULL)
  }
}
n0 <- nchk()
frmtmb_register_frame_check(make_check("refuse_negative_x"))
frmtmb_register_frame_check(make_check("log_only"))
cat("FACTORY checks added:", nchk() - n0, "(2 registered)\n")
set.seed(1)
d <- data.frame(x = rnorm(40), y = rnorm(40))
r <- tryCatch({frm(y ~ x, data = d); "fitted"},
              error = function(e) paste("refused:", conditionMessage(e)))
cat("FACTORY fit with negative x:", r, "\n")
cat("FACTORY checks that ran:", sort(ls(ran)), "\n")

## 2. factory providers ----------------------------------------------------
p0 <- nprv()
make_prov <- function(lb) {
  force(lb)
  function(spec, frame) set_prior("student_t(3, 0, 2.5)", class = "sigma",
                                  lb = lb)
}
frmtmb_register_prior_defaults(make_prov(0.5))
frmtmb_register_prior_defaults(make_prov(0.25))
cat("FACTORY providers added:", nprv() - p0, "(2 registered)\n")
pe <- ns$frmtmb_prior_defaults
cat("FACTORY sigma lb of route = 'sample':",
    tryCatch({
      tab <- default_prior(y ~ x, data = d, route = "sample")
      paste(tab$lb[tab$class == "sigma"], collapse = " ")
    }, error = function(e) conditionMessage(e)), "\n")
ce <- ns$frmtmb_frame_checks
ce$fns <- ce$fns[seq_len(n0)]
pe$providers <- pe$providers[seq_len(p0)]

## 3. library / unload / library -------------------------------------------
gc(); m0 <- sum(gc()[, 2])
t1 <- proc.time()[["elapsed"]]
for (p in c("frmtmb.eam", "frmtmb.ode", "frmtmb.sample")) {
  suppressMessages(library(p, character.only = TRUE))
}
t2 <- proc.time()[["elapsed"]]
gc(); m1 <- sum(gc()[, 2])
cat(sprintf("LOAD eam+ode+sample %.2f s, +%.1f MB; checks=%d providers=%d\n",
            t2 - t1, m1 - m0, nchk(), nprv()))
unloadNamespace("frmtmb.ode")
cat("UNLOAD ode: checks", nchk(), "\n")
suppressMessages(library(frmtmb.ode))
cat("RELOAD ode: checks", nchk(), "\n")
unloadNamespace("frmtmb.sample")
suppressMessages(library(frmtmb.sample))
cat("RELOAD sample: providers", nprv(), "\n")

## 4. cost of a key whose topenv() is frmtmb's namespace -------------------
te <- testthat::test_env("frmtmb")
cat("topenv of a testthat test_env('frmtmb'):",
    environmentName(topenv(te)), "\n")
f <- eval(quote(function(spec, frame) invisible(NULL)), te)
gc(); m0 <- sum(gc()[, 2])
t1 <- proc.time()[["elapsed"]]
frmtmb_register_frame_check(f)
t2 <- proc.time()[["elapsed"]]
gc(); m1 <- sum(gc()[, 2])
cat(sprintf("COST one registration from a test env: %.2f s, +%.1f MB\n",
            t2 - t1, m1 - m0))
t1 <- proc.time()[["elapsed"]]
frmtmb_register_frame_check(f)
t2 <- proc.time()[["elapsed"]]
cat(sprintf("COST again: %.2f s; checks=%d\n", t2 - t1, nchk()))
