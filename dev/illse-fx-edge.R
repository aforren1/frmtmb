# Lane illse fix round 1: degree-6 raw polynomial with (1 | g) on data
# with no group effect (the reviewer's dev/illse-rev-edge4.R design):
# the end state, the warnings, and grad_headroom()'s pieces.
#   Rscript dev/illse-fx-edge.R base|lane|dev deg seed...
args <- commandArgs(TRUE)
arm <- args[1]
R <- "C:/Users/adf44/source/r"
.libPaths(c(if (arm == "lane") file.path(R, "wt-illse-lib"),
            file.path(R, "rellib-r8"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (arm == "dev") {
  suppressMessages(pkgload::load_all(file.path(R, "frmtmb-wt-illse"),
                                     quiet = TRUE))
} else suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
deg <- as.integer(args[2])
for (seed in as.integer(args[-(1:2)])) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = stats::runif(n, 1, 2),
                  g = factor(rep(1:20, length.out = n)))
  d$y <- sin(3 * d$x) + stats::rnorm(n, 0, 0.3)
  xn <- paste0("x", seq_len(deg))
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                                "+ (1 | g)"))
  w <- character()
  f <- withCallingHandlers(frm(fo, data = d), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  lmf <- stats::lm(stats::as.formula(paste("y ~", paste(xn, collapse = "+"))),
                   data = d)
  sdr <- ns$sdr_of(f)
  cat(sprintf("deg %d seed %d | code %d | lm - fit %.2e | polish %s | sd %.1e | root %s | lost %s | explained %s\n",
              deg, seed, f$opt$convergence,
              as.numeric(logLik(lmf)) - as.numeric(logLik(f)),
              paste(signif(f$opt$polish %||% NA, 3), collapse = ","),
              exp(f$opt$par[["theta"]]), !is.null(sdr$cov_fixed_root),
              paste(names(sdr$se_lost), sdr$se_lost, collapse = " "),
              format(f$cache$se_explained %||% "-")))
  cat(paste(" W:", substr(w, 1, 220)), sep = "\n")
  g <- drop(f$obj$gr(f$opt$par))
  cat(" grad:", signif(g, 3), "\n")
  se <- suppressWarnings(sqrt(diag(vcov(f, full = TRUE))))
  cat(" se:", signif(se, 4), "\n")
}
