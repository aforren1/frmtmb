# Reviewer, fix round 1: which data loglik_df() reassembles from, and the
# SE sentence per case. Log dev/catmulti-rev-out/r1b.txt
.libPaths(c("C:/Users/adf44/source/r/wt-catmulti-lib", "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
setwd("C:/Users/adf44/source/r/frmtmb-wt-catmulti")
suppressPackageStartupMessages(library(frmtmb))
q <- function(e) suppressWarnings(suppressMessages(e))
set.seed(5); n <- 400
dh <- data.frame(x = rnorm(n), x2 = rexp(n) + 0.1)
E <- cbind(0, 0.3 + 0.8 * dh$x, -0.4 - 0.5 * log(dh$x2)); P <- exp(E) / rowSums(exp(E))
dh$y <- factor(apply(P, 1, function(p) sample(c("a", "b", "c"), 1, prob = p)))
f <- q(frm(bf(y ~ log(x2)), family = categorical(refcat = NA), data = dh))
cat("names(fit$data):", names(f$data), "\n")
cat("names(frame$data_frame):", names(f$frame$data_frame), "\n")
ll <- logLik(f); cat("log(x2): df", attr(ll, "df"), "why", attr(ll, "df_why"), "\n")
local({ dloc <- dh; g <- q(frm(bf(y ~ x), family = categorical(refcat = NA), data = dloc)); rm(dloc)
  cat("data removed, plain x: df", attr(logLik(g), "df"), "\n") })
x3 <- dh$x
g2 <- q(frm(bf(y ~ x3), family = categorical(refcat = NA), data = dh[, c("y", "x2")]))
cat("variable from the workspace: df", attr(logLik(g2), "df"), "why", attr(logLik(g2), "df_why"), "\n")
fo <- q(frm(bf(y ~ x + offset(x2)), family = categorical(refcat = NA), data = dh))
cat("offset: df", attr(logLik(fo), "df"), "why", attr(logLik(fo), "df_why"), "\n")
fs <- q(frm(bf(y ~ s(x, k = 5)), family = categorical(refcat = NA), data = dh))
cat("s(x): df", attr(logLik(fs), "df"), "why", attr(logLik(fs), "df_why"),
    "| ref:", attr(logLik(q(frm(bf(y ~ s(x, k = 5)), family = categorical(), data = dh))), "df"), "\n")
sentence <- "Under refcat = NA"
w <- function(expr) { out <- character(0); withCallingHandlers(suppressMessages(expr),
  warning = function(cw) { out <<- c(out, conditionMessage(cw)); invokeRestart("muffleWarning") }); out }
cat("dirichlet NA mentions:", any(grepl(sentence, w(frm(bf(Y ~ x), family = dirichlet(refcat = NA), data = local({
  set.seed(2); d <- data.frame(x = rnorm(200)); G <- matrix(rgamma(600, 3), 200)
  d$Y <- G / rowSums(G); colnames(d$Y) <- c("a", "b", "c"); d }))))), "\n")
dd <- data.frame(x = rnorm(100), g = factor(rep(1:10, 10)))
dd$y <- rnorm(100) + 0 * as.integer(dd$g)
wg <- w(frm(bf(y ~ x + (1 | g)), data = dd))
cat("gaussian with sd at 0: warnings", length(wg), "mention:", any(grepl(sentence, wg)), "\n")
dsep <- data.frame(x = c(rnorm(50, -3), rnorm(50, 3)))
dsep$y <- factor(c(rep("a", 50), rep("b", 50)), levels = c("a", "b", "c")); dsep$y[c(1, 60)] <- "c"
dsep$y[dsep$x > 0 & dsep$y == "a"] <- "b"
ws <- w(frm(bf(y ~ x), family = categorical(refcat = "b"), data = dsep))
cat("categorical refcat b separated: warnings", length(ws), "mention:", any(grepl(sentence, ws)), "\n")
cat(substr(ws, 1, 200), sep = "\n")
