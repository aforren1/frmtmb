# How the optimizer fares with link_ndt = "log" against the default
# bounded link. Arm A ("scalar"): one non-decision time, where the two
# links are the same model and must reach the same maximum. Arm B
# ("hier"): ndt ~ 1 + (1 | s), 20 subjects, where brms's log link, the
# default global bound and ndt_group() are three different constraints.
# Usage: Rscript dev/famuni-wiener-ndtrep.R scalar|hier
# Seeds 20261100 + 1:R. Output: dev/famuni-out/wiener-ndtrep-<arm>.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all(
  "C:/Users/adf44/source/r/frmtmb-wt-famuni/extensions/frmtmb.eam",
  quiet = TRUE))
arm <- commandArgs(TRUE)[1]
stopifnot(arm %in% c("scalar", "hier"))

fit_one <- function(form, fam, data) {
  w <- character(0)
  t0 <- proc.time()[["elapsed"]]
  f <- withCallingHandlers(
    tryCatch(frm(form, family = fam, data = data), error = function(e) e),
    warning = function(cw) {
      w <<- c(w, conditionMessage(cw))
      invokeRestart("muffleWarning")
    })
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(f, "error")) {
    return(list(ok = FALSE, err = substr(conditionMessage(f), 1, 80),
                warn = length(w), el = el))
  }
  list(ok = TRUE, fit = f, conv = f$opt$convergence, warn = length(w),
       wtxt = if (length(w)) substr(w[1], 1, 60) else "",
       ll = as.numeric(logLik(f)), el = el)
}

if (arm == "scalar") {
  cat("seed n arm conv warn logLik ndt ndt_se el\n")
  for (n in c(30L, 400L)) {
    for (r in 1:100) {
      seed <- 20261100L + r
      set.seed(seed)
      x <- rnorm(n)
      d <- ddm_simulate(n, mu = 0.8 + 0.5 * x, bs = 1.3, ndt = 0.3,
                        bias = 0.5)
      d$x <- x
      for (lk in c("default", "log")) {
        fam <- if (lk == "default") wiener() else wiener(link_ndt = "log")
        o <- fit_one(bf(rt | dec(upper) ~ x), fam, d)
        if (!o$ok) {
          cat(seed, n, lk, "ERROR", o$warn, NA, NA, NA, o$el, o$err, "\n")
          next
        }
        sp <- summary(o$fit)$spec_pars
        cat(seed, n, lk, o$conv, o$warn, sprintf("%.10f", o$ll),
            sprintf("%.8f", sp["ndt", "Estimate"]),
            sprintf("%.8f", sp["ndt", "Est.Error"]), sprintf("%.2f", o$el),
            o$wtxt, "\n")
      }
    }
  }
} else {
  cat("seed arm conv warn logLik rmse_ms sd_ndt el\n")
  S <- 20L; m <- 50L
  for (r in 1:20) {
    seed <- 20261100L + r
    set.seed(seed)
    s <- factor(rep(seq_len(S), each = m))
    lndt <- log(0.25) + rnorm(S, 0, 0.12)
    ndt_s <- exp(lndt)[as.integer(s)]
    d <- ddm_simulate(S * m, mu = 0.8, bs = 1.4, ndt = ndt_s, bias = 0.5)
    d$s <- s
    arms <- list(
      global = list(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1 + (1 | s),
                       bias = 0.5), wiener()),
      ndt_group = list(bf(rt | dec(upper) + ndt_group(s) ~ 1, bs ~ 1,
                          ndt ~ 1 + (1 | s), bias = 0.5), wiener()),
      log = list(bf(rt | dec(upper) ~ 1, bs ~ 1, ndt ~ 1 + (1 | s),
                    bias = 0.5), wiener(link_ndt = "log")))
    for (nm in names(arms)) {
      o <- fit_one(arms[[nm]][[1]], arms[[nm]][[2]], d)
      if (!o$ok) {
        cat(seed, nm, "ERROR", o$warn, NA, NA, NA, o$el, o$err, "\n")
        next
      }
      nd <- tapply(ndt_time(o$fit), s, mean)
      rmse <- 1000 * sqrt(mean((nd - exp(lndt))^2))
      vc <- tryCatch(frmtmb::VarCorr(o$fit)$s$sd[1, "Estimate"],
                     error = function(e) NA_real_)
      cat(seed, nm, o$conv, o$warn, sprintf("%.6f", o$ll),
          sprintf("%.3f", rmse), sprintf("%.5f", vc), sprintf("%.2f", o$el),
          o$wtxt, "\n")
    }
  }
}
cat("DONE\n")
