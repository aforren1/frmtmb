# Review of lane ci701: D2 beside 0.70.0's other verdicts. A degree-8
# raw polynomial (test-se-flags.R's "illcond" case) with (1 + z | g) and
# no slope variance, so the correlation runs to +/-1; and the rs20 shape
# with check_se = "stop" and "ignore". Lane build, as installed and with
# se_edge_cor() stubbed, in one process. Prints every warning and
# message.
#
#   Rscript dev/ci701-rev-compose.R
.libPaths(c("C:/Users/adf44/source/r/wt-ci701-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
orig <- get("se_edge_cor", ns)
stub <- function(fit, np) logical(np)
conds <- function(expr) {
  out <- character(0)
  v <- withCallingHandlers(tryCatch(expr, error = function(e) {
    out <<- c(out, paste("ERROR:", conditionMessage(e))); NULL
  }), warning = function(w) {
    out <<- c(out, paste("W:", conditionMessage(w)))
    invokeRestart("muffleWarning")
  }, message = function(m) {
    out <<- c(out, paste("M:", conditionMessage(m)))
    invokeRestart("muffleMessage")
  })
  list(value = v, conds = out)
}
show <- function(tag, r) {
  fit <- r$value
  cat("RES", tag, "")
  if (is.null(fit)) {
    cat("no fit\n")
  } else {
    lost <- ns$sdr_of(fit)$se_lost
    th <- fit$estimates$theta
    cat(sprintf("code=%d t=%.4g lost=%s\n", fit$opt$convergence,
                th[length(th)],
                if (length(lost)) paste(names(lost), lost, sep = ":",
                                        collapse = ",") else "none"))
  }
  for (x in r$conds) cat("   ", substr(gsub("\\s+", " ", x), 1, 200), "\n")
}

set.seed(3)
d <- data.frame(x = stats::runif(200, 1, 2), g = factor(rep(1:20, 10)),
                z = stats::rnorm(200))
d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.3)[d$g] +
  stats::rnorm(200, 0, 0.2)
for (k in 1:8) d[[paste0("x", k)]] <- d$x^k
fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + z + (1 + z | g)
for (a in c("lane", "noD2")) {
  assignInNamespace("se_edge_cor", if (a == "lane") orig else stub, ns)
  for (t0 in list(NULL, 500)) {
    st <- if (!is.null(t0)) list(theta = c(0, 0, t0))
    r <- conds({
      m <- frm(fo, data = d, start = st)
      vcov(m)
      ci <- confint_varcorr(m)
      m
    })
    show(sprintf("poly8 %s t0=%s", a, if (is.null(t0)) "def" else t0), r)
  }
}

for (a in c("lane", "noD2")) {
  assignInNamespace("se_edge_cor", if (a == "lane") orig else stub, ns)
  for (s in c(7L, 12L)) {
    set.seed(s)
    g <- factor(rep(seq_len(20), each = 10))
    dd <- data.frame(g = g, x = stats::rnorm(200))
    dd$y <- 1 + 0.5 * dd$x + stats::rnorm(20, 0, 0.7)[g] + stats::rnorm(200)
    for (cs in c("warning", "stop", "ignore")) {
      r <- conds({
        m <- frm(bf(y ~ x + (1 + x | g)), family = gaussian(), data = dd,
                 control = frmtmb_control(check_se = cs))
        vcov(m)
        summary(m)
        ci <- confint_varcorr(m)
        m
      })
      show(sprintf("rs20 s=%d %s check_se=%s", s, a, cs), r)
    }
  }
}
assignInNamespace("se_edge_cor", orig, ns)
