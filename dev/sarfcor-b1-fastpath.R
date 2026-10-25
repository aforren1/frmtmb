# Part B1: frmtmb.sample's derived sd_/cor_ columns for a model with a
# REGISTERED field structure (toyfld, the one test-covstruct-seam.R
# registers) beside a (1 + x | g) block over the SAME grouping factor,
# as car(W, gr = loc) + (1 | loc) is. The old draws_cov_fast() (base
# 18b536fc) kept its own list of field names; the new one asks core's
# frmtmb_covstruct_info(). Both arms run in ONE process, interleaved,
# on the same draws; a control arm runs the new code twice. Seed 21.
# Run: Rscript dev/sarfcor-b1-fastpath.R > dev/sarfcor-out/b1-fastpath.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), "|", find.package("frmtmb.sample"),
    "\n")
ns <- asNamespace("frmtmb.sample")
# the base version of draws_cov_fast(), from the commit this lane is on
old_src <- system2("git", c("-C", "C:/Users/adf44/source/r/frmtmb-wt-sarfcor",
                            "show",
                            "18b536fc:extensions/frmtmb.sample/R/draws-brms.R"),
                   stdout = TRUE)
old_env <- new.env(parent = ns)
exprs <- parse(text = old_src)
for (e in exprs) {
  if (is.call(e) && identical(e[[1L]], as.name("<-")) &&
        identical(e[[2L]], as.name("draws_cov_fast"))) {
    eval(e, old_env)
  }
}
stopifnot(is.function(old_env$draws_cov_fast))

frmtmb_register_covstruct(
  "toyfld",
  lpdf = function(b, theta, blk) {
    K <- exp(2 * theta[1]) * blk[["aux"]][["A"]]
    sum(RTMB::dmvnorm(b, 0, K, log = TRUE))
  },
  vcov = function(theta, blk) {
    matrix(exp(2 * theta[1]), 1, 1,
           dimnames = list("sd(toyfld)", "sd(toyfld)"))
  },
  npar = 1L, sd_idx = 1L,
  natural = function(theta, blk) {
    list(term = "sd(toyfld)", type = "sd", est = theta[1],
         grad = matrix(1, 1, 1))
  },
  draw = function(theta, blk) {
    K <- exp(2 * theta[1]) * blk[["aux"]][["A"]]
    drop(crossprod(chol(K), stats::rnorm(nrow(K))))
  },
  levels_coupled = TRUE,
  term = list(
    args = c("A", "gr"), usage = "toyfld(A, gr = g)",
    parse = function(args, env) list(A = args[["A"]]),
    build = function(spec, gr, lookup) {
      A <- lookup(spec[["A"]], "toyfld()")
      lv <- levels(droplevels(as.factor(gr)))
      list(levels = lv, index = match(as.character(gr), lv),
           aux = list(A = A[lv, lv]))
    }))
set.seed(21)
G <- 30
A <- 0.6^abs(outer(1:G, 1:G, "-"))
dimnames(A) <- list(paste0("L", 1:G), paste0("L", 1:G))
d <- data.frame(g = factor(rep(rownames(A), each = 8), levels = rownames(A)),
                x = rnorm(G * 8))
u <- as.vector(t(chol(0.8 * A)) %*% rnorm(G))
gi <- as.integer(d$g)
bg <- rnorm(G, 0, 0.5)
sg <- rnorm(G, 0, 0.4)
d$y <- 1 + 0.5 * d$x + u[gi] + bg[gi] + sg[gi] * d$x +
  rnorm(nrow(d), 0, 0.6)
fit <- frm(y ~ x + toyfld(A, gr = g) + (1 + x | g), data = d,
           data2 = list(A = A))
s <- suppressMessages(frm_sample(fit, chains = 1, iter = 2000, warmup = 1000,
                                 seed = 21, refresh = 0))
cols_all <- colnames(as.matrix(s))
dv <- ns$draws_cov_cols(fit)
cat("derived names:", dv$names, "\n")
vo <- get("hyp_vals_only", envir = ns)(fit)
ev <- get("hyp_env_vals", envir = ns)(fit, vo$vals, vo$comp)
f_old <- old_env$draws_cov_fast(fit, dv$names, vo, ev)
f_new <- ns$draws_cov_fast(fit, dv$names, vo, ev)
cat("fast map, old:", if (is.null(f_old)) "NULL (every draw through hyp_env_vals())" else
  paste(f_old$names, collapse = " "), "\n")
cat("fast map, new:", if (is.null(f_new)) "NULL" else
  paste(f_new$names, collapse = " "), "\n")
m <- as.matrix(s)
m_in <- m[, setdiff(colnames(m), dv$names), drop = FALSE]
dv_old <- dv
dv_old$fast <- f_old
dv_new <- dv
dv_new$fast <- f_new
v_old <- ns$draws_cov_values(m_in, fit, dv_old)
v_new <- ns$draws_cov_values(m_in, fit, dv_new)
cat("draws:", nrow(m_in), " values identical:", identical(v_old, v_new),
    " max |diff|:", format(max(abs(v_old - v_new))), "\n")
# proc.time() ticks at 10 ms here, so each block repeats the call until
# it runs past 1.2 s, and the time is per call: 2 calls of the old arm,
# 40 of the new and of the control
tm <- function(dvx, k) {
  system.time(for (i in seq_len(k)) {
    ns$draws_cov_values(m_in, fit, dvx)
  })[["elapsed"]] / k
}
res <- data.frame(round = integer(), old = numeric(), new = numeric(),
                  new2 = numeric())
for (r in 1:5) {
  res[r, ] <- c(r, tm(dv_old, 2L), tm(dv_new, 40L), tm(dv_new, 40L))
}
print(res)
cat(sprintf("per call, min over 5 rounds: old %.4f s, new %.4f s, control (new again) %.4f s; old/new %.1f, control ratio %.2f\n",
            min(res$old), min(res$new), min(res$new2),
            min(res$old) / min(res$new), min(res$new2) / min(res$new)))
cat("per draw: old", format(1000 * min(res$old) / nrow(m_in), digits = 3),
    "ms, new", format(1000 * min(res$new) / nrow(m_in), digits = 3), "ms\n")
