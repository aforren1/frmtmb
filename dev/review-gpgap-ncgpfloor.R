# Short check of fix round 2, item 5: ncgp's route read at the floor, as data for the
# merge rule. Not a review of that lane. For each noise: frm() on wt-ncgp-lib, what it
# warns, against the closed form at relative 1e-12 and at brms's absolute 1e-12. r is
# the sampler rule's ratio read in closed form: 1e-12 sdgp^2 over sigma^2 of brms's
# maximum. Seed 8901, n 60.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-ncgpfloor.R ncgpfloor.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")), "\n")
n <- 60
for (A in c(1, 1000)) {
  set.seed(8901)
  x <- round(runif(n, 0, 6), 3)
  e <- rnorm(n)
  u <- sort(unique(x)); J <- match(x, u); D2 <- outer(u, u, "-")^2
  for (s in c(1e-3, 1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 1e-8, 0)) {
    y <- A * (sin(2 * x) + s * e)
    d <- data.frame(x = x, y = y)
    w <- character()
    t0 <- proc.time()[["elapsed"]]
    fr <- withCallingHandlers(
      tryCatch(frm(bf(y ~ gp(x)), family = gaussian(), data = d), error = function(er) {
        cat(sprintf("NF A %g s %.0e ERROR %s\n", A, s, substr(conditionMessage(er), 1, 100)))
        NULL
      }),
      warning = function(wr) {
        w <<- c(w, substr(gsub("[\r\n]+", " ", conditionMessage(wr)), 1, 60))
        invokeRestart("muffleWarning")
      }, message = function(m) invokeRestart("muffleMessage"))
    wall <- proc.time()[["elapsed"]] - t0
    if (is.null(fr)) next
    closed <- function(nug) {
      nll <- function(p) {
        s2 <- exp(2 * p[2])
        K <- s2 * exp(-D2 / (2 * exp(2 * p[3]))); diag(K) <- diag(K) + nug(s2)
        V <- K[J, J]; diag(V) <- diag(V) + exp(2 * p[4])
        R <- tryCatch(chol(V), error = function(e) NULL)
        if (is.null(R)) return(1e10)
        r <- backsolve(R, y - p[1], transpose = TRUE)
        0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * n * log(2 * pi)
      }
      best <- NULL
      for (st in list(c(0, log(A) + 1, 0.3, log(A * max(s, 1e-9))),
                      c(0, log(A) + 1.5, 0.5, log(A * max(s, 1e-9))),
                      c(0, log(A) + 0.7, 0.2, log(A * max(s, 1e-9)) - 1))) {
        o <- nlminb(st, nll, control = list(iter.max = 800, eval.max = 1200))
        if (is.null(best) || o$objective < best$objective) best <- o
      }
      best
    }
    a <- closed(function(s2) 1e-12 * s2)
    b <- closed(function(s2) 1e-12)
    sg <- exp(a$par[4]); sd <- exp(a$par[2])
    r <- 1e-12 * sd^2 / exp(2 * a$par[4])
    cat(sprintf(paste0("NF A %g s %.0e | frm sigma %.3g conv %d ll %.3f (%.1fs) | closed rel: sigma %.3g sdgp %.3g ll %.3f, frm short %.3f | ",
                       "closed brms: sigma %.3g ll %.3f | closed sigma/sdgp %.2e r %.3g thr %.3f | warnings %d: %s\n"),
                A, s, as.numeric(sigma(fr)), fr$opt$convergence, -fr$opt$objective, wall,
                sg, sd, -a$objective, fr$opt$objective - a$objective,
                exp(b$par[4]), -b$objective, sg / sd, r, sqrt(2 / n), length(w),
                paste(unique(w), collapse = " | ")))
  }
}
