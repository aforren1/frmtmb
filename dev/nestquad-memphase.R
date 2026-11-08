# Where the nested tape's memory goes: peak working set after each phase
# of one freeze, in one process.
#   Rscript dev/nestquad-memphase.R <n_outer> <ko> <ki> [gc]
# Seed 20261106, bernoulli, 4 inner levels of 5 rows per outer level.
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-nestquad-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
no <- as.integer(a[1]); ko <- as.integer(a[2]); ki <- as.integer(a[3])
set.seed(20261106)
npo <- 4; rpi <- 5; n <- no * npo * rpi
d <- data.frame(ga = factor(rep(seq_len(no), each = npo * rpi)),
                gb = factor(rep(rep(seq_len(npo), each = rpi), no)),
                x = rnorm(n))
ab <- interaction(d$ga, d$gb, drop = TRUE)
d$y <- rbinom(n, 1, plogis(0.2 + 0.5 * d$x + rnorm(no, 0, 0.7)[d$ga] +
                             rnorm(nlevels(ab), 0, 0.5)[ab]))
fl <- frm(bf(y ~ x + (1 | ga/gb)) + bernoulli(), data = d)
mem <- function(tag) {
  m <- ps::ps_memory_info()
  cat(sprintf("PHASE %-14s rss %7.0f MB peak %7.0f MB\n", tag,
              m[["rss"]] / 2^20, m[["peak_wset"]] / 2^20))
}
gc(); mem("laplace fit")
fr <- fl$frame
nll <- ns$build_objective(fr)
tpl <- fr$par_template
lap <- RTMB::MakeADFun(nll, tpl, random = "b", map = fr$map, silent = TRUE)
plan <- ns$quad_check_blocks(fr)
anc <- ns$quad_nest_anchor(lap, plan, fl$opt$par)
gc(); mem("anchor")
pl <- ns$imp_par_list(tpl, fl$opt$par, fr$map)
qo <- ns$quad_nest_objective(fr, plan, anc, c(ko, ki), pl)
mem("builder")
qo$prime(pl)
mem("prime")
if (length(a) > 3) gc()
otpl <- ns$imp_template(tpl, "b")
t0 <- proc.time()[[3]]
obj <- RTMB::MakeADFun(qo$fn, otpl, silent = TRUE)
cat("tape s", proc.time()[[3]] - t0, "\n")
mem("tape")
sz <- ns$quad_tape_size(obj)
cat("SIZE", sz, "stacked", qo$n_stack, "\n")
t0 <- proc.time()[[3]]
for (i in 1:10) { obj$fn(obj$par); obj$gr(obj$par) }
cat("fn+gr x10 s", proc.time()[[3]] - t0, "\n")
mem("fn+gr")
