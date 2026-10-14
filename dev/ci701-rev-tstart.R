# (ci701-rev-tstart.R: starts t0 = 1000 to 10000) Review of lane ci701: the boundary-correlation fit of test-review-v29.R
# on ulp-moved data (seed 8000 + k, a different draw from the lane's
# 1000 + k), from the default start and from theta = (0, 0, 500).
# On the lane build each fit is made twice in the same process: as
# installed, and with se_edge_cor() stubbed to FALSE (0.70.0's check
# with D1 kept), so that the D2 effect is read without P/E-core noise.
#
#   Rscript dev/ci701-rev-batt.R <lane|base> <k0> <k1>
args <- commandArgs(TRUE)
arm <- args[1]; k0 <- as.integer(args[2]); k1 <- as.integer(args[3])
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
lane <- "C:/Users/adf44/source/r/wt-ci701-lib"
.libPaths(if (arm == "base") c(base, ul) else c(lane, base, ul))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "blas:", extSoftVersion()[["BLAS"]], La_library(), "\n")
set.seed(61)
ng <- 12
u <- stats::rnorm(ng, 0, 0.8)
per <- 5
dd <- data.frame(id = factor(rep(paste0("a", seq_len(ng)), each = per)),
                 x = stats::rnorm(ng * per))
dd$y1 <- 1 + 0.5 * dd$x + u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.4)
dd$y2 <- -1 + 0.2 * dd$x + 0.6 * u[as.integer(dd$id)] +
  stats::rnorm(nrow(dd), 0, 0.3)
d0 <- dd
mdata <- function(k) {
  dd <- d0
  if (k > 0) {
    set.seed(8000 + k)
    for (v in c("x", "y1", "y2")) {
      dd[[v]] <- dd[[v]] * (1 + sample(0:4, nrow(dd), TRUE) * 2^-52)
    }
  }
  dd
}
f <- mvbf(bf(y1 ~ x + (1 | q | id)) + gaussian(),
          bf(y2 ~ x + (1 | q | id)) + gaussian())
orig <- if (exists("se_edge_cor", ns, inherits = FALSE)) get("se_edge_cor", ns)
stub <- function(fit, np) logical(np)
one <- function(k, t0, tag) {
  w <- character(0); bm <- 0L
  st <- if (!is.null(t0)) list(theta = c(0, 0, t0))
  fit <- withCallingHandlers(
    frm(f, data = mdata(k), start = st),
    warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    },
    frmtmb_boundary_fit = function(m) {
      bm <<- bm + 1L; invokeRestart("muffleMessage")
    })
  p <- fit$opt$par
  j <- max(which(names(p) == "theta"))
  q <- p; q[j] <- q[j] + 2 * sign(q[j])
  dl <- fit$obj$fn(p) - fit$obj$fn(q)
  g <- max(abs(drop(fit$obj$gr(p))))
  lost <- ns$sdr_of(fit)$se_lost
  cw <- character(0)
  ci <- withCallingHandlers(confint_varcorr(fit), warning = function(x) {
    cw <<- c(cw, conditionMessage(x)); invokeRestart("muffleWarning")
  })
  s <- ci[ci$type == "sd", ]
  an <- ns$se_fit_analysis(fit, force = TRUE)
  th <- which(names(p) == "theta")
  se <- suppressWarnings(sqrt(diag(an$shown))[th[1:2]])
  cat(sprintf(paste0("RES %s k=%d t0=%s code=%d t=%.1f g=%.3g dl=%.3g ",
                     "nwarn=%d bmsg=%d ciw=%d sdfin=%d lost=%s ",
                     "flost=%s fse=%s w=%s\n"),
              tag, k, if (is.null(t0)) "def" else t0, fit$opt$convergence,
              p[j], g, dl, length(w), bm, length(cw),
              sum(is.finite(c(s$lwr, s$upr))),
              if (length(lost)) paste(names(lost), lost, sep = ":",
                                      collapse = ",") else "none",
              if (length(an$lost)) paste(names(an$lost), an$lost, sep = ":",
                                         collapse = ",") else "none",
              paste(sprintf("%.5f", se), collapse = ","),
              paste(substr(w, 1, 50), collapse = " | ")))
}
for (k in k0:k1) {
  for (t0 in as.list(c(1000, 2000, 3000, 5000, 10000))) {
    if (arm == "base") {
      one(k, t0, "base")
    } else {
      assignInNamespace("se_edge_cor", orig, ns)
      one(k, t0, "lane")
      assignInNamespace("se_edge_cor", stub, ns)
      one(k, t0, "noD2")
      assignInNamespace("se_edge_cor", orig, ns)
    }
  }
}
