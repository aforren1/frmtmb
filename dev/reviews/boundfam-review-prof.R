# Review of lane boundfam: run one test file on the lane build and, at
# every firing of se_edge_limit(), profile each "end" coefficient with the
# OTHER parameters re-optimized (the lane's own probe holds them). A real
# maximum shows as a profile objective that rises again further out.
#
#   Rscript dev/reviews/boundfam-review-prof.R <package> <test file> <out>
#
# Output lines:
#   FIRE <k> <name> est= dir= unit= kinds=
#   PROF <k> <name> t=<natural units along dir> held=<obj - f0> prof=<obj - f0>
local({
  .libPaths(c("C:/Users/adf44/source/r/wt-boundfam-lib",
              "C:/Users/adf44/source/r/rellib-r20",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressMessages(library(frmtmb))
  a <- commandArgs(trailingOnly = TRUE)
  out <- a[3]
  cat("TRACE lib", dirname(find.package("frmtmb")), "\n", file = out)
  ns <- asNamespace("frmtmb")
  rec <- new.env()
  rec$k <- 0L
  assign(".bfrev_rec", rec, envir = globalenv())
  rec$prof <- function(env) {
    if (!exists("lim", envir = env, inherits = FALSE)) return(invisible())
    lim <- get("lim", env)
    if (!length(lim)) return(invisible())
    fit <- get("fit", env); p <- get("p", env); u <- get("u", env)
    dir <- get("dir", env); om <- get("om", env)
    rec$k <- rec$k + 1L
    k <- rec$k
    li <- ns$se_limit_info(fit)
    fam <- tryCatch(fit$spec$responses[[1L]]$family$family,
                    error = function(e) "?")
    obj <- fit$obj
    saved <- ns$obj_state_save(obj)
    on.exit(ns$obj_state_restore(obj, saved), add = TRUE)
    fq <- function(q) {
      v <- tryCatch(suppressWarnings(obj$fn(q)), error = function(e) NA)
      if (is.finite(v)) v else NA_real_
    }
    f0 <- fq(p)
    lw <- fit$opt$lower %||% fit$bounds$lower
    for (j in lim) {
      cat("FIRE ", k, " ", om$names[j], " fam=", fam, " est=",
          signif(p[j], 6), " dir=", dir[j], " unit=", signif(u[j], 4),
          " kinds=", paste(names(li$kind), li$kind, sep = ":",
                           collapse = ","), "\n",
          sep = "", file = out, append = TRUE)
      others <- setdiff(seq_along(p), j)
      for (t in c(-2, 0, 2, 6, 14)) {
        pj <- p[j] + dir[j] * t * u[j]
        q0 <- p; q0[j] <- pj
        held <- fq(q0) - f0
        g <- function(r) {
          q <- p; q[others] <- r; q[j] <- pj
          v <- fq(q)
          if (is.finite(v)) v else 1e100
        }
        gr <- function(r) {
          q <- p; q[others] <- r; q[j] <- pj
          gg <- tryCatch(suppressWarnings(obj$gr(q))[others],
                         error = function(e) rep(NA_real_, length(r)))
          if (all(is.finite(gg))) gg else rep(0, length(r))
        }
        o <- tryCatch(suppressWarnings(stats::nlminb(
          p[others], g, gr, scale = 1 / u[others],
          control = list(iter.max = 400, eval.max = 800))),
          error = function(e) list(objective = NA))
        cat("PROF ", k, " ", om$names[j], " t=", t, " held=",
            signif(held, 6), " prof=", signif(o$objective - f0, 6), "\n",
            sep = "", file = out, append = TRUE)
      }
    }
    invisible()
  }
  suppressMessages(trace("se_edge_limit", where = ns, print = FALSE,
    exit = quote(get(".bfrev_rec", globalenv())$prof(environment()))))
})
Sys.setenv(FRMTMB_RELLIB = "C:/Users/adf44/source/r/wt-boundfam-lib")
source("C:/Users/adf44/source/r/frmtmb-wt-boundfam/dev/release/run-tests.R")
