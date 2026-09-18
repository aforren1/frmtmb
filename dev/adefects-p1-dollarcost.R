# Punch round 1, MAJOR 1. What a `$.frmtmb_fit` method would cost, so
# that D4's choice between "store the frame twice" and "serve it through
# a `$` method" is made on measurements rather than on the assertion the
# lane shipped.
#
#   Rscript dev/adefects-p1-dollarcost.R > dev/adefects-log/p1-dollar.txt 2>&1
#
# PART A counts the `$` reads a fit takes during the operations the
# coordinator named, by REGISTERING a counting `$.frmtmb_fit` and
# running them. The count is only a count if the method is reached from
# inside frmtmb's own namespace as well as from the global environment,
# so that is asserted before anything is counted.
#
# PART B is the per-read cost, on dev/famlink-p2-dollarcost.R's design:
# arms interleaved in ONE process, each block grown past 1.2 s, the
# minimum of several rounds, Sys.time() because proc.time() ticks at
# 10 ms here, and a CONTROL arm built from the same code that must
# report 1.0.
#
# PART C is the on-disk cost of the duplicate, which is what MAJOR 1
# measured: R's serializer does not deduplicate a shared SEXP.
.libPaths(c("C:/Users/adf44/source/r/adefects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "\n")

set.seed(20260917)
n <- 240
d <- data.frame(g = factor(rep(1:12, each = 20)), x = rnorm(n))
d$y <- d$x + rnorm(12, 0, 0.7)[d$g] + rnorm(n)

# ------------------------------------------------------------- PART A

hits <- 0L
`$.frmtmb_fit` <- function(x, name) {
  hits <<- hits + 1L
  # plain `$` on a list partial-matches, so a faithful counter must too
  .subset2(x, name, exact = FALSE)
}
registerS3method("$", "frmtmb_fit", `$.frmtmb_fit`)

# The method is only an instrument if frmtmb's OWN code dispatches to
# it. `fit$frame` inside predict() is the read that decides.
fit0 <- frm(y ~ x + (1 | g), d)
hits <- 0L
invisible(predict(fit0))
reached <- hits > 0L
cat("\n== A. is a registered $.frmtmb_fit reached from inside frmtmb?",
    reached, "(reads seen:", hits, ")\n")
if (!reached) {
  cat("  the count below would be user-level reads only; stopping\n")
  stop("the instrument does not see frmtmb's own reads")
}

count_of <- function(label, expr) {
  hits <<- 0L
  invisible(force(expr))
  cat(sprintf("  %-40s %8d\n", label, hits))
  hits
}
nd <- d[1:20, ]
cat("  `$` reads on a fit, per user action:\n")
ct <- c(
  fit = count_of("frm(y ~ x + (1 | g))", frm(y ~ x + (1 | g), d)),
  predict = count_of("predict(fit)", predict(fit0)),
  predict_nd = count_of("predict(fit, newdata, se.fit = TRUE)",
                        predict(fit0, newdata = nd, se.fit = TRUE)),
  fitted = count_of("fitted(fit)", fitted(fit0)),
  simulate = count_of("simulate(fit, nsim = 5)",
                      simulate(fit0, nsim = 5, seed = 1)),
  summary = count_of("summary(fit)", summary(fit0)),
  ranef = count_of("ranef(fit)", ranef(fit0)),
  vcov = count_of("vcov(fit)", vcov(fit0)),
  boot = count_of("frm_bootstrap(fit, nsim = 20)",
                  frm_bootstrap(fit0, nsim = 20, seed = 2)),
  interop = count_of("insight::get_data / model.frame",
                     model.frame(fit0)))
if (requireNamespace("emmeans", quietly = TRUE)) {
  ct["emmeans"] <- count_of("emmeans::emmeans(fit, 'x')",
                            summary(emmeans::emmeans(fit0, "x")))
}
cat("  the largest single action:", max(ct), "reads\n")

# ------------------------------------------------------------- PART B

# The method this lane would actually ship under option (a2): `data` is
# served from the frame and every other name reads the list.
srv <- function(x, name) {
  if (name == "data") return(.subset2(x, "frame")[["data_frame"]])
  .subset2(x, name, exact = FALSE)
}
`$.adefects_fitcost` <- srv
registerS3method("$", "adefects_fitcost", `$.adefects_fitcost`)
plainlist <- unclass(fit0)
obj <- structure(unclass(fit0), class = c("adefects_fitcost", "list"))
stopifnot(identical(obj$frame, plainlist$frame),
          identical(obj$data, plainlist$frame[["data_frame"]]))

blk <- function(expr) {
  e <- substitute(expr)
  g <- eval(call("function", as.pairlist(alist(n = )),
                 call("for", as.name("i"), quote(seq_len(n)), e)),
            parent.frame())
  function(n) {
    gc(FALSE)
    t0 <- Sys.time()
    g(n)
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
}
arms <- list(
  plain = blk(plainlist$frame),
  plain2 = blk(plainlist$frame),
  method = blk(obj$frame),
  method_data = blk(obj$data),
  bracket = blk(plainlist[["frame"]]))
N <- vapply(arms, function(a) {
  k <- 1e3
  while (a(k) < 1.2) k <- k * 2
  k
}, 1)
ROUNDS <- 5L
best <- setNames(rep(Inf, length(arms)), names(arms))
for (r in seq_len(ROUNDS)) {
  for (a in names(arms)) {
    best[[a]] <- min(best[[a]], arms[[a]](N[[a]]) / N[[a]])
  }
}
cat("\n== B. one `$` read, seconds, minimum of", ROUNDS, "rounds\n")
for (a in names(arms)) {
  cat(sprintf("  %-12s %10.3f us  (block of %g)\n", a,
              best[[a]] * 1e6, N[[a]]))
}
cat(sprintf("  CONTROL plain2 / plain = %.4f (must read about 1.0)\n",
            best[["plain2"]] / best[["plain"]]))
extra <- best[["method"]] - best[["plain"]]
cat(sprintf("  a method read costs %.3f us more than a plain one\n",
            extra * 1e6))
cat(sprintf("  at %d reads, the largest action above, that is %.3f ms\n",
            max(ct), max(ct) * extra * 1e3))

# ------------------------------------------------------------- PART C

cat("\n== C. what the duplicate costs on disk\n")
sizes <- function(obj) {
  raw <- length(serialize(obj, NULL))
  f <- tempfile()
  on.exit(unlink(f), add = TRUE)
  saveRDS(obj, f, compress = "gzip")
  c(raw = raw, gz = file.size(f))
}
for (nn in c(240L, 20000L)) {
  set.seed(20260917)
  dd <- data.frame(g = factor(rep(seq_len(nn / 20), each = 20)),
                   x = rnorm(nn))
  dd$y <- dd$x + rnorm(nn / 20, 0, 0.7)[dd$g] + rnorm(nn)
  f1 <- frm(y ~ x + (1 | g), dd)
  # the fit as this lane ships it, and the same fit with the duplicate
  # element removed, which is what a `$` method would serve instead
  f0 <- f1
  f0$data <- NULL
  s1 <- sizes(f1)
  s0 <- sizes(f0)
  cat(sprintf("  n = %6d  with data: raw %9d gz %8d\n", nn, s1[["raw"]],
              s1[["gz"]]))
  cat(sprintf("               without : raw %9d gz %8d\n", s0[["raw"]],
              s0[["gz"]]))
  cat(sprintf("               delta   : raw %+9d (%+.2f%%) gz %+8d (%+.2f%%)\n",
              s1[["raw"]] - s0[["raw"]],
              100 * (s1[["raw"]] / s0[["raw"]] - 1),
              s1[["gz"]] - s0[["gz"]],
              100 * (s1[["gz"]] / s0[["gz"]] - 1)))
  cat(sprintf("               the frame alone: raw %9d\n",
              length(serialize(f1$frame[["data_frame"]], NULL))))
}
# the control for PART C: a list holding one frame twice, shared, costs
# the same as one holding two equal copies, which is what says the
# serializer does not deduplicate
mf <- data.frame(a = rnorm(20000), b = rnorm(20000))
cp <- data.frame(a = mf$a, b = mf$b)
cat(sprintf("\n  control, shared  list(mf, mf)  raw %9d\n",
            length(serialize(list(mf, mf), NULL))))
cat(sprintf("  control, copies  list(mf, cp)  raw %9d\n",
            length(serialize(list(mf, cp), NULL))))
cat(sprintf("  control, one     list(mf)      raw %9d\n",
            length(serialize(list(mf), NULL))))

cat("\nDONE\n")
