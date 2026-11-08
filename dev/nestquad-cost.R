# Cost of the nested quadrature tape: predicted and actual tape size,
# peak working set and wall time, one configuration per process.
#   Rscript dev/nestquad-cost.R <lib|base> <family> <n_outer>
#           <inner_per_outer> <rows_per_inner> <ko> <ki> [max_gb]
# Prints one COST line. Seed 20261106.
a <- commandArgs(TRUE)
lib <- "C:/Users/adf44/source/r/wt-nestquad-lib"
base <- "C:/Users/adf44/source/r/rellib-r18"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(lib, base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
fam <- a[2]
no <- as.integer(a[3]); npo <- as.integer(a[4]); rpi <- as.integer(a[5])
ko <- as.integer(a[6]); ki <- as.integer(a[7])
max_gb <- if (length(a) >= 8) as.numeric(a[8]) else 4
set.seed(20261106)
n <- no * npo * rpi
d <- data.frame(ga = factor(rep(seq_len(no), each = npo * rpi)),
                gb = factor(rep(rep(seq_len(npo), each = rpi), no)),
                x = rnorm(n))
ab <- interaction(d$ga, d$gb, drop = TRUE)
eta <- 0.2 + 0.5 * d$x + rnorm(no, 0, 0.7)[d$ga] +
  rnorm(nlevels(ab), 0, 0.5)[ab]
d$y <- switch(fam,
  bernoulli = rbinom(n, 1, plogis(eta)),
  poisson = rpois(n, exp(eta)),
  cumulative = factor(cut(eta + rlogis(n), c(-Inf, -0.8, 0.8, Inf),
                          labels = FALSE), ordered = TRUE),
  beta = { p <- plogis(eta); rbeta(n, p * 8, (1 - p) * 8) },
  gaussian = eta + rnorm(n, 0, 0.6))
famobj <- switch(fam, bernoulli = bernoulli(), poisson = poisson(),
                 cumulative = cumulative(), beta = Beta(),
                 gaussian = gaussian())
peak0 <- ps::ps_memory_info()[["peak_wset"]]
msgs <- character(0)
t0 <- proc.time()[[3]]
f <- tryCatch(withCallingHandlers(
  frm(bf(y ~ x + (1 | ga/gb)), family = famobj, data = d,
      quadrature = TRUE, verbose = TRUE,
      control = frmtmb_control(quad_nodes = c(ko, ki),
                               quad_max_gb = max_gb)),
  message = function(m) {
    msgs <<- c(msgs, conditionMessage(m))
    invokeRestart("muffleMessage")
  }, warning = function(w) invokeRestart("muffleWarning")),
  error = function(e) e)
el <- proc.time()[[3]] - t0
peak1 <- ps::ps_memory_info()[["peak_wset"]]
pred <- regmatches(msgs, regexpr("about [0-9.e+]+ MB predicted", msgs))
tape_s <- grep("nested quadrature tape", msgs, value = TRUE)
sz <- if (!inherits(f, "error")) frmtmb:::quad_tape_size(f$obj)
tape_mb <- if (is.null(sz)) NA else
  sum(sz * frmtmb:::quad_nest_entry_bytes) / 2^20
cat(sprintf(paste("COST fam %s outer %d per %d rows %d n %d nodes %dx%d",
                  "time %.2fs peak0MB %.0f peakMB %.0f deltaMB %.0f",
                  "pred %s tapeMB %.1f ops %s values %s inputs %s %s\n"),
            fam, no, npo, rpi, n, ko, ki, el, peak0 / 2^20, peak1 / 2^20,
            (peak1 - peak0) / 2^20,
            if (length(pred)) sub("about ([0-9.e+]+) MB.*", "\\1", pred[1])
            else NA, tape_mb,
            if (is.null(sz)) NA else sz[1], if (is.null(sz)) NA else sz[2],
            if (is.null(sz)) NA else sz[3],
            if (inherits(f, "error")) paste("ERROR", conditionMessage(f))
            else sprintf("logLik %.8f", as.numeric(logLik(f)))))
writeLines(msgs)
