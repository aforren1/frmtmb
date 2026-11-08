# Where an exgaussian nested quadrature fit at 3000 rows ends the R
# process: the anchor, each probe tape, or the guard.
source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/review-nestquad-common.R")
n <- as.integer(commandArgs(TRUE)[1] %||% "3000")
set.seed(20261108)
n_out <- n %/% 20
d <- data.frame(school = factor(rep(seq_len(n_out), each = 20)),
                cls = rep(rep(1:4, each = 5), n_out), x = rnorm(n))
d$cls_uni <- factor(paste0(d$school, "_", d$cls))
eta <- 0.3 + 0.4 * d$x + rnorm(n_out, 0, 0.5)[d$school] +
  rnorm(nlevels(d$cls_uni), 0, 0.5)[d$cls_uni]
d$y <- eta + rnorm(n) + rexp(n, 2)
fl <- frm(bf(y ~ x + (1 | school/cls_uni)) + exgaussian(), data = d)
fr <- fl$frame
say <- function(...) { cat(sprintf("%s peak %.0f cur %.0f\n", paste0(...), peak_mb(), cur_mb())); flush(stdout()) }
say("laplace done")
lap <- RTMB::MakeADFun(ns$build_objective(fr), fr$par_template, random = "b",
                       map = fr$map, silent = TRUE)
plan <- ns$quad_check_blocks(fr); say("plan")
anc <- ns$quad_nest_anchor(lap, plan, fl$opt$par); say("anchor")
pl <- ns$imp_par_list(fr$par_template, fl$opt$par, fr$map)
otpl <- ns$imp_template(fr$par_template, "b")
for (k in list(c(1L, 1L), c(1L, 2L), c(2L, 1L))) {
  qo <- ns$quad_nest_objective(fr, plan, anc, k, pl)
  say("objective ", k[1], "x", k[2], " n_stack ", qo$n_stack)
  obj <- RTMB::MakeADFun(qo$fn, otpl, silent = TRUE)
  say("taped ", k[1], "x", k[2])
  print(ns$quad_tape_size(obj))
  rm(obj); gc()
}
b <- ns$quad_nest_guard(fr, plan, anc, c(11L, 15L), pl, Inf, otpl, NULL)
say("guard ", format(b / 2^20), " MB")
cat("DONE\n")
