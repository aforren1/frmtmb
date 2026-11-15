# Short check of parity2's fix round, F3: brms 2.23.0's default display
# names on the shapes of parity2-short-a.R. Arm "internal": brm(empty =
# TRUE) and brms's own get_all_effects() + prepare_cond_data(reorder =
# TRUE), the steps conditional_effects.brmsfit() takes to name its
# displays (no compile). Arm "fit": fixed_param fits, the real
# conditional_effects() names.
#   Rscript dev/reviews/parity2-short-brms.R internal|fit > <out> 2>&1
arm <- commandArgs(TRUE)[1]
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
options(mc.cores = 1)
cat("brms", as.character(packageVersion("brms")), "arm", arm, "\n")
set.seed(21)
n <- 200
e <- data.frame(x = rnorm(n), z = rnorm(n),
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                h = factor(sample(c("u", "v"), n, TRUE)),
                grp = factor(sample(letters[1:8], n, TRUE)),
                kk = sample(1:6, n, TRUE))
e$ch <- as.character(e$f)
e$lg <- e$h == "u"
e$y <- 1 + e$x + 0.5 * e$z + as.integer(e$f) * 0.3 + rnorm(8)[e$grp] +
  rnorm(6)[e$kk] + rnorm(n)
shapes <- list(
  fx = y ~ f * x, xf = y ~ x * f, fh = y ~ f * h, xz = y ~ x * z,
  fx_re = y ~ f * x + (1 | grp), sby = y ~ s(x, by = f),
  chx = y ~ ch * x, lgx = y ~ lg * x, kkx_grp = y ~ kk * x + (1 | kk))
named <- function(nm) {
  if (nm %in% c("fx", "xf")) return(list("f:x", "x:f"))
  if (nm == "chx") return(list("ch:x", "x:ch"))
  if (nm == "lgx") return(list("lg:x", "x:lg"))
  list()
}
internal_names <- function(fit, effects = NULL) {
  bt <- brms:::brmsterms(fit$formula)
  rsv <- brms:::rsv_vars(bt)
  def <- is.null(effects)
  eff <- if (def) brms:::get_all_effects(bt, rsv_vars = rsv) else
    strsplit(effects, ":")
  mf <- model.frame(fit)
  cond <- data.frame(row.names = 1)
  vapply(eff, function(ev) {
    cd <- brms:::prepare_cond_data(
      mf[, ev, drop = FALSE], conditions = cond,
      int_vars = brms:::get_int_vars(bt),
      group_vars = brms:::get_group_vars(bt), resolution = 5,
      reorder = def)
    paste(attr(cd, "effects"), collapse = ":")
  }, "")
}
for (nm in names(shapes)) {
  if (identical(arm, "internal")) {
    fit <- brm(shapes[[nm]], data = e, empty = TRUE)
    cat(nm, ":", internal_names(fit), "\n")
    for (ef in named(nm)) cat("  effects =", ef, ":",
                              internal_names(fit, ef), "\n")
  } else {
    fit <- suppressMessages(suppressWarnings(brm(
      shapes[[nm]], data = e, algorithm = "fixed_param", chains = 1,
      iter = 20, warmup = 0, init = 0, refresh = 0, seed = 1, silent = 2)))
    ce <- suppressWarnings(conditional_effects(fit, resolution = 5))
    cat(nm, ":", names(ce), "\n")
    for (ef in named(nm)) {
      cat("  effects =", ef, ":",
          names(suppressWarnings(conditional_effects(fit, effects = ef,
                                                     resolution = 5))),
          "\n")
    }
    k <- grep(":", names(ce), value = TRUE)
    if (length(k)) cat("  ", k, "effect1__ class",
                       class(ce[[k]]$effect1__)[1], "\n")
  }
}
cat("done\n")
