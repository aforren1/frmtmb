# brms 2.23.0 fits for lane catmulti: categorical(refcat =), multinomial()
# bare and with refcat, and the simplex families inside mvbf(). Records
# the post-fit shapes and dimnames, and saves each fit for the
# agreement checks. Output: dev/catmulti-out/brms-fits.txt and
# dev/catmulti-out/brms-<name>.rds (local).
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
out <- "dev/catmulti-out"
source("dev/catmulti-data.R")
sec <- function(t) cat("\n========", t, "========\n")
shape <- function(nm, x) {
  cat(nm, ": dim", paste(dim(x), collapse = " x "), "\n")
  dn <- dimnames(x)
  for (k in seq_along(dn)) {
    cat("  dimnames[[", k, "]]: ",
        if (is.null(dn[[k]])) "NULL" else
          paste(utils::head(dn[[k]], 6), collapse = ", "), "\n", sep = "")
  }
}
post <- function(fit, nm) {
  saveRDS(fit, file.path(out, paste0("brms-", nm, ".rds")))
  print(fixef(fit))
  cat("max Rhat", max(rhat(fit), na.rm = TRUE), "\n")
  for (f in c("fitted", "posterior_epred", "posterior_linpred",
              "predict", "posterior_predict", "log_lik")) {
    r <- tryCatch(do.call(f, list(fit)),
                  error = function(e) paste("ERROR:", conditionMessage(e)))
    if (is.character(r)) cat(f, ":", r, "\n") else shape(f, r)
  }
  r <- tryCatch(fitted(fit, scale = "linear"),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r)) cat("fitted linear :", r, "\n") else
    shape("fitted(scale = 'linear')", r)
  r <- tryCatch(posterior_linpred(fit, transform = TRUE),
                error = function(e) paste("ERROR:", conditionMessage(e)))
  if (is.character(r)) cat("linpred transform :", r, "\n") else
    shape("posterior_linpred(transform = TRUE)", r)
  invisible(fit)
}
ctl <- list(chains = 2, iter = 2000, warmup = 1000, refresh = 0,
            backend = "rstan", seed = 1)
run <- function(nm, ...) {
  sec(nm)
  fit <- tryCatch(do.call(brm, c(list(...), ctl)),
                  error = function(e) {
                    cat("ERROR:", conditionMessage(e), "\n"); NULL
                  })
  if (!is.null(fit)) post(fit, nm)
}
dc <- catmulti_cat_data(11)
run("cat-refb", bf(f ~ x), data = dc, family = categorical(refcat = "b"))
run("cat-refna", bf(f ~ x), data = dc, family = categorical(refcat = NA))
fc1 <- readRDS(file.path(out, "brms-cat-refb.rds"))
sec("cat-refb posterior_linpred head")
pl <- posterior_linpred(fc1)
print(pl[1:2, 1:3, ])
print(fitted(fc1)[1:2, , ])
dm <- catmulti_mult_data(12)
run("mult-bare", bf(Yc | trials(N) ~ x), data = dm, family = multinomial())
run("mult-refy2", bf(Yc | trials(N) ~ x), data = dm,
    family = multinomial(refcat = "y2"))
fm <- readRDS(file.path(out, "brms-mult-bare.rds"))
sec("mult-bare fitted head")
print(fitted(fm)[1:2, , ])
print(posterior_linpred(fm)[1:2, 1:3, ])
ds <- catmulti_mv_data(13)
run("mv-dir", bf(Y ~ x, family = dirichlet()) +
      bf(w ~ x, family = gaussian()) + set_rescor(FALSE), data = ds)
run("mv-ln", bf(Y ~ x, family = logistic_normal()) +
      bf(w ~ x, family = gaussian()) + set_rescor(FALSE), data = ds)
