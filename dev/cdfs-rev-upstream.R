# Reviewer: reproduce RTMB-8 (pbeta shape derivative at whole shapes),
# RTMB-9 (pgamma third derivatives; upper tail not taped), and check
# whether log_ibeta_half()'s pbeta() region (a b / (a + b) > 150) meets
# RTMB-8 through a whole SECOND shape, which is how negbinomial() calls
# it (b = q + 1). brms-32 and brms-33 are in cdfs-rev-brms31.R.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("RTMB", format(packageVersion("RTMB")), "\n")
fd_pb <- function(x, a, b, j) {
  h <- 1e-6
  if (j == 1) (pbeta(x, a + h, b) - pbeta(x, a - h, b)) / (2 * h) else
    (pbeta(x, a, b + h) - pbeta(x, a, b - h)) / (2 * h)
}
cat("== RTMB-8: d pbeta(x, a, b) / d(a, b)\n")
for (cs in list(c(0.3, 2, 5), c(0.3, 2.0001, 5), c(0.3, 2.5, 5), c(0.3, 2.5, 5.5),
                c(0.3, 2, 5.5), c(0.4, 200, 300), c(0.4, 200.5, 300),
                c(0.4, 200.5, 300.5), c(0.4, 200, 300.5))) {
  x <- cs[1]
  tp <- RTMB::MakeTape(function(p) RTMB::pbeta(x, p[1], p[2]), cs[2:3])
  ad <- as.numeric(tp$jacobian(cs[2:3]))
  cat(sprintf("x %.1f a %-8g b %-8g: AD %s  FD %s\n", x, cs[2], cs[3],
              paste(signif(ad, 6), collapse = ", "),
              paste(signif(c(fd_pb(x, cs[2], cs[3], 1), fd_pb(x, cs[2], cs[3], 2)), 6),
                    collapse = ", ")))
}
cat("== log_ibeta_half() in its pbeta() region, whole b (negbinomial's q + 1)\n")
for (cs in list(c(0.6, 300.3, 201), c(0.6, 300, 201), c(0.5, 500.7, 500))) {
  x <- cs[1]
  tp <- RTMB::MakeTape(function(p) ns$log_ibeta_half(x, p[1], cs[3]), cs[2])
  ad <- as.numeric(tp$jacobian(cs[2]))
  fd <- (pbeta(x, cs[2] + 1e-6, cs[3], log.p = TRUE) -
           pbeta(x, cs[2] - 1e-6, cs[3], log.p = TRUE)) / 2e-6
  cat(sprintf("x %.2f a %-7g b %-5g: d/da AD %.8g FD %.8g\n", x, cs[2], cs[3], ad, fd))
}
cat("== negbinomial lcdf: d/d log shape at a whole shape, large shape and count\n")
for (cs in list(c(250, 300, 260), c(250, 300.5, 260), c(1000, 400, 1010), c(5, 2, 4))) {
  q <- cs[1]
  f <- function(p) ns$nbinom_cdf_log(q, log(cs[3]), p[1])$lF
  ad <- as.numeric(RTMB::MakeTape(f, log(cs[2]))$jacobian(log(cs[2])))
  h <- 1e-6
  fd <- (pnbinom(q, size = cs[2] * exp(h), mu = cs[3], log.p = TRUE) -
           pnbinom(q, size = cs[2] * exp(-h), mu = cs[3], log.p = TRUE)) / (2 * h)
  cat(sprintf("q %g shape %-6g mu %g: AD %.8g FD %.8g\n", q, cs[2], cs[3], ad, fd))
}
cat("== RTMB-9: pgamma third derivatives; the upper tail on the tape\n")
set.seed(1)
T3 <- RTMB::MakeTape(function(p) log(RTMB::pgamma(exp(p[1]), shape = exp(p[2]))),
                     c(0, 0))$jacfun()$jacfun()
pts <- cbind(runif(500, -3, 4), runif(500, -2, 6))
nf <- sum(apply(pts, 1, function(p) any(!is.finite(T3(p)))))
cat("non-finite third derivatives at", nf, "of 500 points\n")
r <- tryCatch(RTMB::MakeTape(function(p) RTMB::pgamma(2, shape = p[1], lower.tail = FALSE), 3),
              error = function(e) conditionMessage(e))
cat("pgamma(lower.tail = FALSE) on the tape:", if (is.character(r)) r else "tapes", "\n")
r <- tryCatch(RTMB::MakeTape(function(p) RTMB::ppois(2, p[1], lower.tail = FALSE), 3),
              error = function(e) conditionMessage(e))
cat("ppois(lower.tail = FALSE) on the tape:", if (is.character(r)) r else "tapes", "\n")
