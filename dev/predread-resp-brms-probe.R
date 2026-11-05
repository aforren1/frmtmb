# brms 2.23.0 on every reading of `resp` the lane touches, RUN (data seed 5).
#   Rscript dev/predread-resp-brms-probe.R > dev/predread-log/resp-brms-probe.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
quiet <- function(e) suppressWarnings(suppressMessages(e))
tr <- function(lab, e) {
  cat("\n---", lab, "\n")
  r <- tryCatch(quiet(e), error = function(err) {
    cat("ERROR:", conditionMessage(err), "\n")
    invisible(NULL)
  })
  if (!is.null(r)) {
    if (is.list(r) && !is.data.frame(r)) {
      cat("names:", names(r), "\n")
      cat("class:", class(r), "\n")
      if (length(r)) cat("first names:", names(r[[1L]]), "\n")
    } else if (!is.null(dim(r))) {
      cat("dim:", dim(r), "\n")
      if (!is.null(dimnames(r))) {
        cat("rownames:", utils::head(rownames(r), 4), "\n")
        if (length(dim(r)) == 3L) cat("third:", dimnames(r)[[3L]], "\n")
      }
    } else {
      print(r)
    }
  }
  invisible(r)
}
set.seed(5)
n <- 80
d <- data.frame(x = rnorm(n))
d$y <- d$x + rnorm(n)
d$y2 <- -d$x + rnorm(n)
d$y3 <- 0.5 * d$x + rnorm(n)
u <- quiet(brm(bf(y ~ x), data = d, chains = 1, iter = 300, refresh = 0,
               seed = 1))
for (a in list(list(lab = "character(0)", r = character(0)),
               list(lab = "c('y', 'y')", r = c("y", "y")),
               list(lab = "'y'", r = "y"),
               list(lab = "'nosuch'", r = "nosuch"),
               list(lab = "c('nosuch', 'other')", r = c("nosuch", "other")),
               list(lab = "NA_character_", r = NA_character_),
               list(lab = "1L", r = 1L))) {
  tr(paste("uni posterior_epred resp =", a$lab),
     posterior_epred(u, resp = a$r))
  tr(paste("uni posterior_linpred resp =", a$lab),
     posterior_linpred(u, resp = a$r))
  tr(paste("uni posterior_predict resp =", a$lab),
     posterior_predict(u, resp = a$r))
  tr(paste("uni fitted resp =", a$lab), fitted(u, resp = a$r))
  tr(paste("uni predict resp =", a$lab), predict(u, resp = a$r))
  tr(paste("uni log_lik resp =", a$lab), log_lik(u, resp = a$r))
  tr(paste("uni bayes_R2 resp =", a$lab), bayes_R2(u, resp = a$r))
  tr(paste("uni conditional_effects resp =", a$lab),
     conditional_effects(u, "x", resp = a$r))
}

mv <- quiet(brm(mvbf(bf(y ~ x), bf(y2 ~ x), bf(y3 ~ x), rescor = FALSE),
                data = d, chains = 1, iter = 300, refresh = 0, seed = 1))
tr("mv conditional_effects, no resp", conditional_effects(mv, "x"))
tr("mv conditional_effects, no effects", conditional_effects(mv))
tr("mv conditional_effects resp = 'y2'",
   conditional_effects(mv, "x", resp = "y2"))
tr("mv conditional_effects resp = c('y', 'y3')",
   conditional_effects(mv, "x", resp = c("y", "y3")))
tr("mv conditional_effects resp = c('y3', 'y')",
   conditional_effects(mv, "x", resp = c("y3", "y")))
tr("mv conditional_effects resp = c('y', 'y')",
   conditional_effects(mv, "x", resp = c("y", "y")))
tr("mv conditional_effects resp = character(0)",
   conditional_effects(mv, "x", resp = character(0)))
tr("mv conditional_effects resp = 'nosuch'",
   conditional_effects(mv, "x", resp = "nosuch"))
tr("mv conditional_effects, method posterior_predict",
   conditional_effects(mv, "x", method = "posterior_predict"))
tr("mv conditional_effects, several, method posterior_predict",
   conditional_effects(mv, "x", resp = c("y", "y3"),
                       method = "posterior_predict"))
ce <- conditional_effects(mv, "x", resp = c("y", "y3"))
cat("attrs of a panel:", names(attributes(ce[[1L]])), "\n")
cat("columns:", names(ce[[1L]]), "\n")
cat("attr response:", attr(ce[[1L]], "response"), "\n")
cat("attr effects:", attr(ce[[1L]], "effects"), "\n")
tr("mv bayes_R2 resp = character(0)", bayes_R2(mv, resp = character(0)))
tr("mv bayes_R2 resp = c('y', 'y')", bayes_R2(mv, resp = c("y", "y")))
tr("mv loo_R2 resp = character(0)", loo_R2(mv, resp = character(0)))
tr("mv loo_R2 resp = c('y', 'y')", loo_R2(mv, resp = c("y", "y")))
tr("mv posterior_epred resp = c('y', 'y')",
   posterior_epred(mv, resp = c("y", "y")))
tr("mv posterior_epred resp = character(0)",
   posterior_epred(mv, resp = character(0)))

# a multivariate mixture: pp_mixture reads resp
mm <- quiet(brm(mvbf(bf(y ~ x), bf(y2 ~ x), rescor = FALSE),
                data = d, family = mixture(gaussian(), gaussian()),
                chains = 1, iter = 200, refresh = 0, seed = 1, init = 0))
tr("mv mixture pp_mixture no resp", pp_mixture(mm))
tr("mv mixture pp_mixture resp = 'y2'", pp_mixture(mm, resp = "y2"))
tr("mv mixture pp_mixture resp = character(0)",
   pp_mixture(mm, resp = character(0)))
tr("mv mixture pp_mixture resp = c('y', 'y')",
   pp_mixture(mm, resp = c("y", "y")))
tr("mv mixture pp_mixture resp = c('y', 'y2')",
   pp_mixture(mm, resp = c("y", "y2")))
tr("mv mixture pp_mixture resp = 'nosuch'",
   pp_mixture(mm, resp = "nosuch"))
um <- quiet(brm(bf(y ~ x), data = d, family = mixture(gaussian(), gaussian()),
                chains = 1, iter = 200, refresh = 0, seed = 1, init = 0))
tr("uni mixture pp_mixture resp = character(0)",
   pp_mixture(um, resp = character(0)))
tr("uni mixture pp_mixture resp = c('y', 'y')",
   pp_mixture(um, resp = c("y", "y")))
tr("uni mixture pp_mixture resp = 'nosuch'", pp_mixture(um, resp = "nosuch"))
# what a mixed mv with cs() names its layers
d$o <- cut(d$x + rlogis(n), c(-Inf, -1, 1, Inf), labels = FALSE)
ms <- quiet(brm(mvbf(bf(y ~ x), bf(o ~ x + cs(y3), family = sratio()),
                     rescor = FALSE),
                data = d, chains = 1, iter = 200, refresh = 0, seed = 1))
nd <- data.frame(x = c(-1, 1), y3 = 0.2, y = 0, o = 1)
tr("mv + cs posterior_linpred", posterior_linpred(ms, newdata = nd))
tr("mv + cs fitted linear", fitted(ms, newdata = nd, scale = "linear"))
tr("mv + cs emmeans", as.data.frame(emmeans::emmeans(
  ms, ~ x | rep.meas, at = list(x = c(-1, 1)))))
tr("mv + cs posterior_linpred resp = 'o'",
   posterior_linpred(ms, newdata = nd, resp = "o"))
# hurdle_cumulative with cs()
set.seed(6)
dh <- data.frame(x = rnorm(200), z = rnorm(200))
dh$y <- pmax(0L, cut(dh$x + 0.5 * dh$z + rlogis(200),
                     c(-Inf, -1, 0, 1.5, Inf), labels = FALSE) - 1L)
hc <- quiet(brm(bf(y ~ x + cs(z)), data = dh, family = hurdle_cumulative(),
                chains = 1, iter = 200, refresh = 0, seed = 1, init = 0))
tr("hurdle_cumulative cs posterior_linpred",
   posterior_linpred(hc, newdata = data.frame(x = c(-1, 1), z = 0)))
