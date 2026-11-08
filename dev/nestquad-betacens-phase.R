# Peak working set by phase for the Beta + cens fuzz spec, 5 x 5 nodes.
#   Rscript dev/nestquad-betacens-phase.R
.libPaths(c("C:/Users/adf44/source/r/wt-nestquad-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
ns <- asNamespace("frmtmb")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-nestquad/tests/testthat"
env <- new.env(parent = ns)
for (h in c("helper-warnings.R", "helper-fuzz.R")) {
  sys.source(file.path(wt, h), envir = env)
}
plan <- env$fuzz_plan(seed = 20260901L, size = 300L)
sp <- plan[plan$mode == "quadrature" & plan$re == "nested" &
             plan$family == "Beta", , drop = FALSE][1, ]
d <- env$fuzz_data(sp)
dd <- if (is.data.frame(d)) d else d$data
bform <- eval(parse(text = env$fuzz_bf_text(sp)))
fl <- suppressWarnings(frm(bform, data = dd))
mem <- function(tag) {
  m <- ps::ps_memory_info()
  cat(sprintf("PHASE %-10s rss %7.0f MB peak %7.0f MB\n", tag,
              m[["rss"]] / 2^20, m[["peak_wset"]] / 2^20))
}
gc(); mem("laplace")
fr <- fl$frame
tpl <- fr$par_template
lap <- RTMB::MakeADFun(ns$build_objective(fr), tpl, random = "b",
                       map = fr$map, silent = TRUE)
pq <- ns$quad_check_blocks(fr)
anc <- ns$quad_nest_anchor(lap, pq, fl$opt$par)
pl <- ns$imp_par_list(tpl, fl$opt$par, fr$map)
gc(); mem("anchor")
qo <- ns$quad_nest_objective(fr, pq, anc, c(5, 5), pl)
qo$prime(pl)
mem("prime")
gc(); mem("gc")
t0 <- proc.time()[[3]]
obj <- RTMB::MakeADFun(qo$fn, ns$imp_template(tpl, "b"), silent = TRUE)
cat("makeadfun s", proc.time()[[3]] - t0, "\n")
mem("makeadfun")
cat("SIZE", ns$quad_tape_size(obj), "\n")
