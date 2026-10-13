# Reviewer of the 0.70.0 consolidation: lane moparity's new mo()
# structures under lane optrap's units and Newton steps and lane
# seflags' SE check, which no lane saw together. Fits the same designs
# on the release build (r8) and on lane moparity's build (lane), seeds
# 1..N, the response at scale 1 and 1e3, and writes one row per fit.
#   Rscript dev/relrev070-interplay.R r8|lane N out.tsv
args <- commandArgs(TRUE)
arm <- args[1]
N <- as.integer(args[2])
out <- args[3]
.libPaths(c(switch(arm, r8 = "C:/Users/adf44/source/r/rellib-r8",
                   lane = "C:/Users/adf44/source/r/wt-moparity-lib"),
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("arm", arm, as.character(packageVersion("frmtmb")),
    find.package("frmtmb"), "\n")
mk <- function(seed, n = 300, n_g = 20) {
  set.seed(seed)
  dd <- data.frame(
    x = sample(0:3, n, TRUE),
    f = factor(sample(c("a", "b", "c"), n, TRUE)),
    g = factor(sample(paste0("g", seq_len(n_g)), n, TRUE)),
    z = stats::rnorm(n), w = stats::rnorm(n))
  eff <- c(0, 1, 1.5, 1.7)
  bg <- stats::rnorm(n_g, 0, 0.5)
  ag <- stats::rnorm(n_g, 0, 0.4)
  dd$y <- 1 + ag[dd$g] +
    (1 + (dd$f == "b") - (dd$f == "c") + bg[dd$g]) * eff[dd$x + 1] +
    stats::rnorm(n, 0, 0.5)
  dd$cnt <- stats::rpois(n, exp(0.2 + 0.3 * eff[dd$x + 1] *
                                  (1 + (dd$f == "b"))))
  # x never takes the value 2: the simplex's middle step is not read
  # by the data, the case a Dirichlet prior identifies
  dd$x2 <- ifelse(dd$x == 2, 3L, dd$x)
  dd
}
designs <- list(
  mof = list(f = y ~ mo(x) * f, fam = gaussian()),
  mof_g = list(f = y ~ mo(x) * f + (1 | g), fam = gaussian()),
  mo_re = list(f = y ~ mo(x) + (mo(x) | g), fam = gaussian()),
  mozw = list(f = y ~ mo(x) * z * w, fam = gaussian()),
  mof_pois = list(f = cnt ~ mo(x) * f, fam = poisson(), noscale = TRUE),
  dir_map = list(f = y ~ mo(x2), fam = gaussian(),
                 prior = set_prior("dirichlet(c(2, 2, 2))", class = "simo",
                                   coef = "mox21")),
  dir_flat = list(f = y ~ mo(x2), fam = gaussian())
)
rows <- list()
for (seed in seq_len(N)) {
  d0 <- mk(seed)
  for (nm in names(designs)) {
    ds <- designs[[nm]]
    for (sc in if (isTRUE(ds$noscale)) 1 else c(1, 1e3)) {
      d <- d0
      d$y <- d$y * sc
      w <- character()
      fit <- withCallingHandlers(
        tryCatch(do.call(frm, c(list(bf(ds$f), family = ds$fam, data = d),
                                if (!is.null(ds$prior))
                                  list(prior = ds$prior))),
                 error = function(e) e),
        warning = function(x) {
          w <<- c(w, conditionMessage(x))
          invokeRestart("muffleWarning")
        }, message = function(x) invokeRestart("muffleMessage"))
      if (inherits(fit, "error")) {
        rows[[length(rows) + 1L]] <- data.frame(
          seed = seed, design = nm, scale = sc, ll = NA, code = NA,
          polish = NA, mo_runs = NA, kappa = NA, spread_mo = NA,
          nse = NA, nse_nan = NA, nwarn = NA,
          warn = substr(conditionMessage(fit), 1, 200))
        next
      }
      fe <- suppressWarnings(fixef(fit))
      se <- fe[, 2]
      lp <- fit$frame$linpreds[[1L]]
      sp <- if (exists("predictor_col_spread", ns)) {
        s <- get("predictor_col_spread", ns)(lp)
        mc <- vapply(lp[["mo"]] %||% list(), function(m) as.integer(m$col), 0L)
        paste(s[mc], collapse = ",")
      } else NA
      kap <- if (exists("predictor_design_kappa", ns)) {
        get("predictor_design_kappa", ns)(fit$frame)
      } else NA
      # summary() runs the SE check; collect its warnings too
      w2 <- character()
      withCallingHandlers(summary(fit), warning = function(x) {
        w2 <<- c(w2, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
      ww <- unique(c(w, w2))
      rows[[length(rows) + 1L]] <- data.frame(
        seed = seed, design = nm, scale = sc,
        ll = as.numeric(logLik(fit)) + nrow(d) * log(sc) *
          (!isTRUE(ds$noscale)),
        code = fit$opt$convergence,
        polish = if (is.null(fit$opt$polish)) 0 else fit$opt$polish[["steps"]],
        mo_runs = if (is.null(fit$opt$mo_search)) 0 else
          fit$opt$mo_search[["runs"]],
        kappa = kap, spread_mo = sp, nse = length(se),
        nse_nan = sum(!is.finite(se)), nwarn = length(ww),
        warn = substr(gsub("[\t\n]", " ", paste(ww, collapse = " | ")), 1,
                      300))
    }
  }
}
X <- do.call(rbind, rows)
utils::write.table(X, out, sep = "\t", quote = FALSE, row.names = FALSE)
cat("fits", nrow(X), "\n")
