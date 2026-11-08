# Brief step 6 (bigx x nestquad): a shifted covariate under nested
# quadrature fits as the covariate centered. Nested rounds run in natural
# coordinates; the final tape carries the design's basis.
#
#   Rscript dev/rel080-nestshift.R <lib> [seeds]
args <- commandArgs(trailingOnly = TRUE)
lib <- args[1L]
seeds <- if (length(args) > 1L) eval(parse(text = args[2L])) else 1:3
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")), "\n")
sim <- function(seed, n_s = 15, n_c = 3, n_r = 4) {
  set.seed(seed)
  d <- expand.grid(r = seq_len(n_r), class = seq_len(n_c), school = seq_len(n_s))
  d$school <- factor(d$school)
  d$class <- factor(paste(d$school, d$class))
  us <- rnorm(n_s, 0, 0.8)[d$school]
  uc <- rnorm(nlevels(d$class), 0, 0.8)[d$class]
  d$x <- rnorm(nrow(d))
  d$y <- rbinom(nrow(d), 1, plogis(0.2 + 0.7 * d$x + us + uc))
  d
}
for (s in seeds) {
  d <- sim(s)
  for (m in c(0, 1e4, 1e6)) {
    d$w <- d$x + m
    t0 <- proc.time()[["elapsed"]]
    w <- character()
    f <- withCallingHandlers(
      frm(y ~ w + (1 | school/class), data = d, family = bernoulli(),
          quadrature = TRUE),
      warning = function(cnd) {
        w <<- c(w, paste0(class(cnd)[1L], ": ", substr(conditionMessage(cnd), 1, 60)))
        invokeRestart("muffleWarning")
      })
    el <- proc.time()[["elapsed"]] - t0
    sm <- fixef(f)
    cat(sprintf("NS seed %d m %-6g logLik %.8f slope %.7f se %.7f conv %d basis %s nested %s %.1fs warn %d %s\n",
                s, m, as.numeric(logLik(f)), sm["w", 1], sm["w", 2],
                f$opt$convergence,
                !is.null(f$obj$env$.frmtmb_opt_basis), !is.null(f$quad_nested),
                el, length(w), paste(w, collapse = " | ")))
  }
}
cat("NESTSHIFT DONE\n")
