# Fix round, MAJOR 1: the reviewer's tight-tolerance sweep, uninstrumented.
# Designs and the exact reference come from gen() and ref_ll() of
# dev/reviews/quadseam-review-acc.R (read, not copied). Each fit at each
# grad_tol: the reported logLik minus the exact one at the fit's own
# estimate, the settle record, every warning class, seconds.
#   Rscript dev/quadseam-fr-sweep.R <lib|base> <design,...> <seeds> <tols>
# log dev/quadseam-log/fr-sweep-<tag>.txt
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
designs <- strsplit(a[2], ",")[[1]]
seeds <- eval(parse(text = a[3]))
tols <- as.numeric(strsplit(a[4], ",")[[1]])
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
src <- readLines(
  "C:/Users/adf44/source/r/frmtmb-wt-quadseam/dev/reviews/quadseam-review-acc.R")
eval(parse(text = src[grep("^gen <- function", src):
                        (grep("^for \\(des in designs\\)", src) - 1L)]))
for (des in designs) for (s in seeds) for (tol in tols) {
  g <- gen(des, s)
  fam <- g$fam
  w <- character()
  t0 <- proc.time()[[3]]
  f <- tryCatch(withCallingHandlers(
    frm(fo[[fam]], family = frm_fam[[fam]], data = g$d, quadrature = TRUE,
        control = frmtmb_control(grad_tol = tol)),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    },
    message = function(x) invokeRestart("muffleMessage")),
    error = function(e) e)
  tt <- proc.time()[[3]] - t0
  if (inherits(f, "error")) {
    cat(sprintf("FIT %s %d tol %g ERROR %s\n", des, s, tol,
                substr(conditionMessage(f), 1, 100)))
    next
  }
  beta <- unname(f$estimates$beta); sd <- exp(unname(f$estimates$theta[1]))
  shape <- if (fam == "gam") exp(unname(f$estimates$betad[1])) else NA
  ex <- ref_ll(fam, g$d, beta, sd, shape)
  qs <- f$opt[["quad_settle"]]
  cat(sprintf(paste0("FIT %s %d tol %g ll %.10f err %10.3e rounds %s ",
                     "last %s grad %s warn %s secs %.1f\n"),
              des, s, tol, as.numeric(logLik(f)),
              as.numeric(logLik(f)) - ex,
              if (is.null(qs)) "-" else qs[["rounds"]],
              if (is.null(qs) || is.na(qs["last"])) "-" else
                format(qs[["last"]], digits = 3),
              if (is.null(qs)) "-" else format(qs[["grad"]], digits = 3),
              if (length(w)) paste(w, collapse = ",") else "-", tt))
}
cat("END\n")
