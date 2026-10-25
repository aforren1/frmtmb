# Reviewer's brms 2.23.0 probes of bhaz(), no compile.
#   Rscript dev/review-bhaz-brms.R > dev/review-bhaz-log/brms.txt 2>&1
.libPaths(c("C:/Users/adf44/source/r/wt-bhaz-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb.survival)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.survival"),
    "brms", format(packageVersion("brms")), "splines2",
    format(packageVersion("splines2")), "\n")
cat("brms::cox formals:", names(formals(brms::cox)), "\n")
cat("splines2::mSpline formals:", names(formals(splines2::mSpline)), "\n")
set.seed(11); n <- 200; x <- rnorm(n)
tt <- rexp(n, exp(-0.5 + 0.7 * x)); ct <- rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x,
                 g = factor(sample(c("a", "b", "c"), n, TRUE)))
cat("range time", format(range(dd$time), digits = 10), "\n")
try1 <- function(lab, expr) {
  w <- character()
  r <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(w2) {
      w <<- c(w, conditionMessage(w2)); invokeRestart("muffleWarning")
    })
  if (inherits(r, "error")) {
    cat(sprintf("%-40s ERROR: %s\n", lab, conditionMessage(r)))
  } else {
    cat(sprintf("%-40s OK\n", lab))
  }
  if (length(w)) cat("   warnings:", unique(w), sep = "\n     ")
  if (!inherits(r, "error")) r else NULL
}
sdat <- function(fo, d = dd) {
  brms::standata(brms::bf(fo), data = d, family = brms::cox())
}
show <- function(sd) {
  if (is.null(sd)) return(invisible())
  cat(sprintf("   Kbhaz %d  knots %s  boundary %s  ngrbhaz %s\n",
              sd$Kbhaz, paste(format(attr(sd$Zbhaz, "knots"), digits = 8),
                              collapse = " "),
              paste(format(attr(sd$Zbhaz, "Boundary.knots"), digits = 8),
                    collapse = " "),
              if (is.null(sd$ngrbhaz)) "-" else sd$ngrbhaz))
}
r <- range(dd$time); d <- diff(r)
cat("claimed default boundary c(max(min - d/50, 0), max + d/50):",
    format(c(max(r[1] - d / 50, 0), r[2] + d / 50), digits = 10), "\n")
show(try1("bhaz() default", sdat(time | cens(cens) ~ x)))
show(try1("bhaz(df = 6)", sdat(time | cens(cens) + bhaz(df = 6) ~ x)))
show(try1("bhaz(knots, Boundary.knots)",
          sdat(time | cens(cens) + bhaz(knots = c(0.5, 1.1),
                                        Boundary.knots = c(0, 16)) ~ x)))
show(try1("bhaz(df = 8, knots = 2 knots)",
          sdat(time | cens(cens) + bhaz(df = 8, knots = c(0.5, 1.1)) ~ x)))
show(try1("bhaz(gr = g)", sdat(time | cens(cens) + bhaz(gr = g) ~ x)))
sg <- try1("bhaz(g) unnamed", sdat(time | cens(cens) + bhaz(g) ~ x))
show(sg)
show(try1("bhaz(6) unnamed", sdat(time | cens(cens) + bhaz(6) ~ x)))
show(try1("bhaz(df = 4.5)", sdat(time | cens(cens) + bhaz(df = 4.5) ~ x)))
show(try1("bhaz(df = 3)", sdat(time | cens(cens) + bhaz(df = 3) ~ x)))
show(try1("bhaz(df=3, intercept=FALSE)",
          sdat(time | cens(cens) + bhaz(df = 3, intercept = FALSE) ~ x)))
show(try1("bhaz(foo = 1)", sdat(time | cens(cens) + bhaz(foo = 1) ~ x)))
show(try1("bhaz(periodic = TRUE)",
          sdat(time | cens(cens) + bhaz(periodic = TRUE) ~ x)))
show(try1("bhaz(Boundary = c(0, 16)) partial",
          sdat(time | cens(cens) + bhaz(Boundary = c(0, 16)) ~ x)))
show(try1("bhaz(d = 6) partial",
          sdat(time | cens(cens) + bhaz(d = 6) ~ x)))
show(try1("bhaz(deg = 2) partial",
          sdat(time | cens(cens) + bhaz(deg = 2) ~ x)))
try1("brms::cox(bhaz = list(df = 6))", brms::cox(bhaz = list(df = 6)))
try1("bhaz() with gaussian",
     brms::standata(brms::bf(time | bhaz(df = 6) ~ x), data = dd,
                    family = gaussian()))

# Boundary knots inside the data: how bad is the extrapolation?
neg <- function(bk) {
  sd <- try1(paste0("Boundary.knots = c(", paste(bk, collapse = ", "), ")"),
             sdat(time | cens(cens) + bhaz(Boundary.knots = !!bk) ~ x))
  if (is.null(sd)) return(invisible())
  Z <- sd$Zbhaz; C <- sd$Zcbhaz
  out <- dd$time < bk[1] | dd$time > bk[2]
  cat(sprintf(paste0("   rows outside %d; rows of Zbhaz with a negative ",
                     "entry %d; rows with a negative Zbhaz %%*%% (1/K) %d;",
                     " min over rows of max_k Zbhaz %.4g;",
                     " Zcbhaz negative entries %d\n"),
              sum(out), sum(apply(Z < 0, 1, any)),
              sum(Z %*% rep(1 / ncol(Z), ncol(Z)) < 0),
              min(apply(Z, 1, max)), sum(C < 0)))
}
# brms evaluates the flag in the formula's environment; build the
# formulas with the value inlined
for (bk in list(c(0.5, 10), c(0, 14), c(0, 14.4), c(0.01, 14.71))) {
  fo <- eval(bquote(time | cens(cens) + bhaz(Boundary.knots = .(bk)) ~ x))
  sd <- try1(paste0("Boundary.knots = c(", paste(bk, collapse = ", "), ")"),
             sdat(fo))
  if (is.null(sd)) next
  Z <- sd$Zbhaz; C <- sd$Zcbhaz
  out <- dd$time < bk[1] | dd$time > bk[2]
  cat(sprintf(paste0("   rows outside %d; rows of Zbhaz with a negative ",
                     "entry %d; rows with a negative uniform-weight hazard ",
                     "%d; Zcbhaz negative entries %d\n"),
              sum(out), sum(apply(Z < 0, 1, any)),
              sum(Z %*% rep(1 / ncol(Z), ncol(Z)) < 0), sum(C < 0)))
}

# frmtmb's basis against brms's standata() for the four spellings, from
# frmtmb's own frame
strip <- function(m) matrix(as.numeric(m), nrow(m))
cmp <- function(lab, fo) {
  sd <- sdat(fo)
  fr <- frm(frmtmb::bf(fo), family = frmtmb.survival::cox(), data = dd,
            dry_run = "frame")
  av <- fr$aterm_values$time
  cat(sprintf("%-32s K %d/%d  M %.3e  I %.3e  bnd brms %s frm %s\n", lab,
              sd$Kbhaz, ncol(av$Zbhaz),
              max(abs(av$Zbhaz - strip(sd$Zbhaz))) / max(abs(sd$Zbhaz)),
              max(abs(av$Zcbhaz - strip(sd$Zcbhaz))) / max(abs(sd$Zcbhaz)),
              paste(format(attr(sd$Zbhaz, "Boundary.knots"), digits = 8),
                    collapse = ","),
              paste(format(av$bhaz_spec$boundary, digits = 8),
                    collapse = ",")))
}
cmp("default", time | cens(cens) ~ x)
cmp("bhaz(df = 7)", time | cens(cens) + bhaz(df = 7) ~ x)
cmp("bhaz(knots = 3, no boundary)",
    time | cens(cens) + bhaz(knots = c(0.3, 1, 3)) ~ x)
cmp("bhaz(degree = 1, df = 4)",
    time | cens(cens) + bhaz(degree = 1, df = 4) ~ x)
cmp("bhaz(Boundary.knots = c(0, 30))",
    time | cens(cens) + bhaz(Boundary.knots = c(0, 30)) ~ x)
