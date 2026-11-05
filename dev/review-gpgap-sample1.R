# Review of lane gpgap, items 2 and 5: frm_sample() against brms, drawn.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-sample.R <log> low <s>
#     low noise: does the sampled field carry the floor on sigma, silently?
#   bash dev/review-gpgap-run.sh dev/review-gpgap-sample.R <log> order
#     ordinary noise on unsorted x with ties: row-wise epred, new positions,
#     draw names. Anything that assumed sorted x would show here.
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
a <- commandArgs(trailingOnly = TRUE)
arm <- a[1]
set.seed(6601)
n <- 40
x <- if (arm == "low") round(runif(n, 0, 6), 3) else
  sample(seq(0, 6, by = 0.2), n, TRUE)
e0 <- rnorm(n)
s <- if (arm == "low") as.numeric(a[2]) else 0.3
d <- data.frame(x = x, y = sin(2 * x) + s * e0)
iter <- as.integer(if (length(a) >= 3) a[3] else 1500L)
warned <- character()
tm <- function(expr) {
  t0 <- proc.time()[["elapsed"]]
  v <- withCallingHandlers(expr, warning = function(w) {
    warned <<- c(warned, paste0(class(w)[1], ": ",
                                substr(conditionMessage(w), 1, 110)))
    invokeRestart("muffleWarning")
  })
  cat(sprintf("  wall %.0f s\n", proc.time()[["elapsed"]] - t0))
  v
}
cat("frm_sample formula route\n")
fs <- tm(frm_sample(bf(y ~ gp(x)), data = d, family = gaussian(),
                    chains = 2, iter = iter, warmup = iter %/% 2, seed = 66,
                    refresh = 0, cores = 1,
                    control = list(adapt_delta = 0.95, max_treedepth = 12)))
cat("warnings from frm_sample:\n")
print(unique(warned))
nug <- any(grepl("frmtmb_gp_nugget", warned))
warned <- character()
cat("brms\n")
bm <- tm(brms::brm(brms::bf(y ~ gp(x)), data = d, family = gaussian(),
                   chains = 2, iter = iter, warmup = iter %/% 2, seed = 66,
                   refresh = 0, cores = 1, silent = 2,
                   control = list(adapt_delta = 0.95, max_treedepth = 12)))
cat("warnings from brm:\n")
print(unique(warned))
ds <- posterior::as_draws_df(fs)
db <- posterior::as_draws_df(bm)
cat("frm_sample draw names:\n")
print(utils::head(posterior::variables(ds), 12))
cat("brms draw names:\n")
print(utils::head(posterior::variables(db), 12))
q <- function(v) signif(stats::quantile(v, c(0.05, 0.5, 0.95)), 4)
rh <- function(dd, nm) signif(posterior::rhat(posterior::extract_variable_matrix(dd, nm)), 4)
sgn <- intersect(c("sigma"), posterior::variables(ds))
cat(sprintf("frm_sample has floor? nugget warning %s\n", nug))
for (pr in list(c("sigma", "sigma"),
                c(grep("^sdgp", posterior::variables(ds), value = TRUE)[1],
                  grep("^sdgp", posterior::variables(db), value = TRUE)[1]),
                c(grep("^lscale", posterior::variables(ds), value = TRUE)[1],
                  grep("^lscale", posterior::variables(db), value = TRUE)[1]))) {
  cat(sprintf("%-22s frm_sample q05/50/95 %s rhat %s | brms (%s) %s rhat %s\n",
              pr[1], paste(q(ds[[pr[1]]]), collapse = " "), rh(ds, pr[1]), pr[2],
              paste(q(db[[pr[2]]]), collapse = " "), rh(db, pr[2])))
}
cat(sprintf("true sigma %.2e; ratio of medians, sigma: %.4g\n", s,
            stats::median(ds$sigma) / stats::median(db$sigma)))
cat("the floor would be sigma^2 + 1e-6 sdgp^2 = total: frm_sample median of ",
    "sqrt(sigma^2 + 1e-6 * sdgp^2) ",
    signif(stats::median(sqrt(ds$sigma^2 + 1e-6 *
      ds[[grep("^sdgp", posterior::variables(ds), value = TRUE)[1]]]^2)), 4),
    "\n")
if (arm == "order") {
  es <- brms::posterior_epred(fs)
  eb <- brms::posterior_epred(bm)
  zr <- (colMeans(es) - colMeans(eb)) / sqrt(apply(eb, 2, var))
  cat(sprintf("row-wise epred at training rows: max |mean diff|/post sd %.3f\n",
              max(abs(zr))))
  cat("corr of frm_sample's epred means with y in data order:",
      signif(cor(colMeans(es), d$y), 4), "\n")
  nd <- data.frame(x = c(0.33, 2.71, 5.55, 6.4))
  es2 <- brms::posterior_epred(fs, newdata = nd)
  eb2 <- brms::posterior_epred(bm, newdata = nd)
  print(signif(rbind(frm = colMeans(es2), brms = colMeans(eb2),
                     brms_sd = apply(eb2, 2, sd)), 4))
  # the first b[] column is the smallest x, not the first row
  b1 <- grep("^b\\[", posterior::variables(ds), value = TRUE)
  if (length(b1)) {
    cat("columns b[...]:", length(b1), "; distinct x:", length(unique(d$x)), "\n")
    mb <- colMeans(as.matrix(ds[, b1]))
    u <- sort(unique(d$x))
    cat("sorted reading, cor(b, sin(2 u)):", signif(cor(mb, sin(2 * u)), 4),
        "; first-appearance reading:",
        signif(cor(mb, sin(2 * unique(d$x))), 4), "\n")
  }
}
