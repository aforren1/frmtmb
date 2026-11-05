# Profile, second pass, per build: where one fit's time goes (Rprof with
# memory profiling), one evaluation's cost (a perturbed point, since a
# repeated point is served from RTMB's cache), tape operators where RTMB
# reports them, and a design with many groups, where the proposal's own
# arithmetic (per group and per outer parameter) is largest.
#   Rscript dev/impsamp-prof2.R <label> [lib]
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
      f(k)
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
  t0 <- proc.time()[["elapsed"]]
  io <- ns$build_importance_objective(fr, lay, gm, plan)
  otpl <- ns$imp_template(tpl, "b")
  obj <- RTMB::MakeADFun(io$fn, otpl, silent = TRUE)
  t_tape <- proc.time()[["elapsed"]] - t0
  peak <- sum(gc()[, 6L])
  info <- tryCatch(TMB:::info(obj$env$ADFun), error = function(e) NULL)
  ops <- tryCatch(
    .Call("InfoADFunObject", obj$env$ADFun$ptr, PACKAGE = "RTMB"),
    error = function(e) NULL)
  jit <- function(k) at + (k %% 2L) * 1e-7
  tf <- clock(function(k) obj$fn(jit(k)))
  tg <- clock(function(k) obj$gr(jit(k)))
  cat(sprintf(
    "%s %s N=%d plan %.3f tape %.3f peak %.0f Mb fn %.5f gr %.5f ops %s\n",
    label, lab, N, t_plan, t_tape, peak, tf, tg,
    if (is.null(ops)) "NA" else paste(unlist(ops)[1:min(4, length(unlist(ops)))],
                                      collapse = "/")))
}
db <- E$eq_dat("binomial", seed = 124)
dp <- E$eq_dat("poisson", seed = 69)
set.seed(99)
G <- 1000L
g <- factor(rep(seq_len(G), each = 5L))
x <- rnorm(5L * G)
dm <- data.frame(y = rpois(5L * G, exp(0.5 + 0.3 * x + rnorm(G, 0, 0.7)[g] +
                                         rnorm(G, 0, 0.3)[g] * x)),
                 x = x, g = g)
for (r in 1:2) {
  cat(sprintf("%s control %.3f s\n", label, control()))
  for (N in c(200L, 1000L)) {
    one("binomial-int", bf(y | trials(ntr) ~ x + (1 | g)), binomial(), db, N)
    one("poisson-slope", bf(y ~ x + (1 + x | g)), poisson(), dp, N)
  }
  one("poisson-slope-1000g", bf(y ~ x + (1 + x | g)), poisson(), dm, 64L)
}
# one whole fit, profiled
tf <- tempfile()
Rprof(tf, memory.profiling = TRUE, interval = 0.01)
fit <- suppressWarnings(frm(bf(y ~ x + (1 + x | g)), family = poisson(),
                            data = dp, importance = 512L))
Rprof(NULL)
s <- summaryRprof(tf, memory = "both")
cat("total", s$sampling.time, "s; rounds", fit$importance$rounds, "\n")
top <- utils::head(s$by.total[order(-s$by.total$total.time), ], 25)
print(top[, c("total.time", "total.pct", "self.time", "mem.total")])
cat("nlminb evaluations in the last round:", fit$opt$evaluations, "\n")
