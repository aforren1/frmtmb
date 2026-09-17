## Reviewer check for lane wt-famlink: what does the exact `$` method
## cost? Counts calls to `$.frmtmb_family` during representative work
## (load-independent), and times one call against `[[` on the same
## object, interleaved, minimum of 7 rounds of 2e5 calls each, with a
## control (`[[` against `[[`) that must report about 1.0.
## Seed 20260916. Usage: Rscript dev/famlink-rev-dollarcost.R
ARM <- "lane"
source("dev/famlink-rev-common.R")
set.seed(20260916)
d <- data.frame(x = rnorm(300), g = factor(rep(1:30, 10)))
d$y <- rpois(300, exp(0.5 + 0.3 * d$x + rnorm(30, 0, 0.3)[d$g]))
ns <- asNamespace("frmtmb")
orig <- get("$.frmtmb_family", ns)
cnt <- new.env(); cnt$n <- 0L
counting <- function(x, name) { cnt$n <- cnt$n + 1L; orig(x, name) }
unlockBinding("$.frmtmb_family", ns); assign("$.frmtmb_family", counting, ns)
registerS3method("$", "frmtmb_family", counting, envir = ns)
tm <- function(label, expr) {
  cnt$n <- 0L
  t <- system.time(expr)[["elapsed"]]
  cat(sprintf("%-40s $ calls %7d  elapsed %.2f s\n", label, cnt$n, t))
}
tm("control: 3 direct f$lpdf reads", { ff <- negbinomial(); ff$lpdf; ff$links; ff$link })
tm("frm negbinomial + (1|g)", fit <- frm(y ~ x + (1 | g), data = d, family = negbinomial()))
tm("predict(newdata)", predict(fit, newdata = d))
tm("simulate(nsim = 50)", simulate(fit, nsim = 50, seed = 1))
tm("frm_bootstrap(nsim = 20)", frm_bootstrap(fit, nsim = 20, seed = 1))
assign("$.frmtmb_family", orig, ns)
registerS3method("$", "frmtmb_family", orig, envir = ns)
f <- negbinomial()
N <- 2e5
one <- function(which) {
  t0 <- proc.time()[["elapsed"]]
  if (which == "dollar") for (i in seq_len(N)) f$lpdf else
    for (i in seq_len(N)) f[["lpdf"]]
  proc.time()[["elapsed"]] - t0
}
r <- replicate(7, c(dollar = one("dollar"), b1 = one("bracket"), b2 = one("bracket")))
cat(sprintf("per call: $ %.2f us, [[ %.2f us; ratio $/[[ %.2f; control [[/[[ %.2f\n",
            1e6 * min(r["dollar", ]) / N, 1e6 * min(r["b1", ]) / N,
            min(r["dollar", ]) / min(r["b1", ]), min(r["b2", ]) / min(r["b1", ])))
