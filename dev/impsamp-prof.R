# Profile the hot loop: one evaluation of the corrected objective and its
# gradient, per build. Load-independent counts first (tape operators and
# the peak of R's heap while it is built), then the clock, with a control
# built from fixed arithmetic. Interleave the arms by running this file
# once per build at the same time:
#   Rscript dev/impsamp-prof.R <label> [lib]
args <- commandArgs(TRUE)
label <- args[[1]]
lib <- if (length(args) > 1L) args[[2]] else ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
E <- new.env(parent = ns)
sys.source("dev/equiv-specs.R", envir = E)
cat("lib:", find.package("frmtmb"), label, "\n")
control <- function() {
  t0 <- proc.time()[["elapsed"]]
  s <- 0
  for (i in 1:3e7) s <- s + i %% 7
  proc.time()[["elapsed"]] - t0
}
clock <- function(f) {
  best <- Inf
  for (r in 1:3) {
    k <- 0L
    t0 <- proc.time()[["elapsed"]]
    repeat {
      f()
      k <- k + 1L
      el <- proc.time()[["elapsed"]] - t0
      if (el > 1.2) break
    }
    best <- min(best, el / k)
  }
  best
}
one <- function(lab, form, fam, d, N) {
  fit <- suppressWarnings(suppressMessages(frm(form, family = fam, data = d)))
  fr <- fit$frame
  tpl <- fit$obj$env$parList(fit$opt$par)
  nll <- ns$build_objective(fr)
  lap <- RTMB::MakeADFun(nll, tpl, random = "b", map = fr[["map"]],
                         silent = TRUE)
  lay <- ns$imp_layout(fr[["re_blocks"]])
  gm <- ns$imp_group_map(fr, lay)
  at <- fit$opt$par
  invisible(gc(reset = TRUE))
  t0 <- proc.time()[["elapsed"]]
  plan <- ns$imp_plan(lap, fr, lay, at, N, 1L)
  t_plan <- proc.time()[["elapsed"]] - t0
  io <- ns$build_importance_objective(fr, lay, gm, plan)
  otpl <- ns$imp_template(tpl, "b")
  t0 <- proc.time()[["elapsed"]]
  obj <- RTMB::MakeADFun(io$fn, otpl, silent = TRUE)
  t_tape <- proc.time()[["elapsed"]] - t0
  peak <- sum(gc()[, 6L])
  ops <- tryCatch({
    tp <- RTMB::MakeTape(function(p) {
      pl <- ns$imp_par_list(otpl, stats::setNames(p, names(at)))
      "[<-" <- RTMB::ADoverload("[<-")
      io$fn(pl)
    }, unname(at))
    v <- utils::capture.output(print(tp))
    paste(v, collapse = " ")
  }, error = function(e) paste("tape count failed:", conditionMessage(e)))
  tf <- clock(function() obj$fn(at))
  tg <- clock(function() obj$gr(at))
  cat(sprintf(
    "%s %s N=%d plan %.3f s tape %.3f s Rheap_peak %.0f Mb fn %.5f s gr %.5f s\n  %s\n",
    label, lab, N, t_plan, t_tape, peak, tf, tg, substr(ops, 1, 160)))
}
db <- E$eq_dat("binomial", seed = 124)
dp <- E$eq_dat("poisson", seed = 69)
for (r in 1:2) {
  cat(sprintf("%s control %.3f s\n", label, control()))
  for (N in c(200L, 1000L)) {
    one("binomial-int", bf(y | trials(ntr) ~ x + (1 | g)), binomial(), db, N)
    one("poisson-slope", bf(y ~ x + (1 + x | g)), poisson(), dp, N)
  }
}
