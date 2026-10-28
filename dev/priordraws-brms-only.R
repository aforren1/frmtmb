# brms 2.23.0 on the rest of the questions: the ar prior draw name,
# sample_prior = "only" (what it refuses, what prior_draws() and
# posterior_predict() give), save_warmup = FALSE with inc_warmup = TRUE,
# hypothesis() with and without prior draws, update(sample_prior = ).
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(brms))
set.seed(1)
d <- data.frame(x = rnorm(60), z = rnorm(60), g = factor(rep(1:6, 10)),
                h = factor(rep(1:5, each = 12)), t = rep(1:12, 5))
d$y <- 1 + 0.5 * d$x + rnorm(60)
q <- function(lab, expr) {
  cat("\n=====", lab, "=====\n")
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n"); NULL
    }),
    warning = function(w) {
      cat("WARNING:", conditionMessage(w), "\n")
      invokeRestart("muffleWarning")
    })
  invisible(r)
}
fitq <- function(...) suppressMessages(brm(..., data = d, chains = 1,
                                            iter = 200, seed = 1,
                                            refresh = 0))
q("AR stancode", {
  sc <- stancode(y ~ x + ar(time = t, gr = h), data = d,
                 prior = prior(normal(0, 0.5), class = ar),
                 sample_prior = "yes")
  cat(grep("ar|prior_", strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
})
q("AR default only", stancode(y ~ x + ar(time = t, gr = h), data = d,
                              prior = prior(normal(0, 1), class = b),
                              sample_prior = "only"))
q("smooth only", stancode(y ~ s(x), data = d,
                          prior = prior(normal(0, 1), class = b),
                          sample_prior = "only"))
q("sd flat only", stancode(y ~ x + (1 | g), data = d,
                           prior = c(prior(normal(0, 1), class = b),
                                     set_prior("", class = "sd")),
                           sample_prior = "only"))
q("only stancode likelihood", {
  sc <- stancode(y ~ x + (1 | g), data = d,
                 prior = prior(normal(0, 1), class = b),
                 sample_prior = "only")
  cat(grep("prior_only|target", strsplit(sc, "\n")[[1]], value = TRUE),
      sep = "\n")
})
fo <- q("only fit", fitq(y ~ x + (1 | g),
                         prior = prior(normal(0, 1), class = b),
                         sample_prior = "only"))
if (!is.null(fo)) {
  q("only prior_draws", print(dim(prior_draws(fo))))
  q("only variables", cat(variables(fo), "\n"))
  q("only pp", print(dim(posterior_predict(fo))))
  q("only summary b_x", print(posterior_summary(fo, variable = "b_x")))
  q("only hypothesis", print(hypothesis(fo, "x = 0")$hypothesis))
}
fy <- q("yes fit", fitq(y ~ x + (1 | g),
                        prior = prior(normal(0, 1), class = b),
                        sample_prior = "yes", save_warmup = FALSE))
if (!is.null(fy)) {
  q("nowarm sim", print(fy$fit@sim[c("iter", "warmup", "warmup2",
                                     "n_save")]))
  q("nowarm inc_warmup df", print(nrow(as_draws_df(fy, inc_warmup = TRUE))))
  q("nowarm as.mcmc inc_warmup",
    print(dim(as.mcmc(fy, inc_warmup = TRUE)[[1]])))
  q("nowarm nsamples", print(nsamples(fy, incl_warmup = TRUE)))
  q("hyp x", print(hypothesis(fy, "x = 0", seed = 1)$hypothesis))
  q("hyp Intercept", print(hypothesis(fy, "Intercept = 0")$hypothesis))
  q("hyp sigma", print(hypothesis(fy, "sigma = 1", class = NULL)$hypothesis))
  q("hyp sd", print(hypothesis(fy, "sd_g__Intercept = 1",
                               class = NULL)$hypothesis))
  q("hyp two", print(hypothesis(fy, "x - Intercept = 0")$hypothesis))
  q("prior_draws b_x head", print(head(prior_draws(fy, "b_x"), 3)))
  q("prior_draws sigma", print(names(prior_draws(fy, c("sigma", "b_x",
                                                     "sd_g__Intercept")))))
  q("prior_draws Intercept", print(names(prior_draws(fy, "b_Intercept"))))
}
fn <- q("no-prior fit", fitq(y ~ x, prior = prior(normal(0, 1), class = b)))
if (!is.null(fn)) {
  q("no-prior prior_draws", print(class(prior_draws(fn))))
  q("no-prior prior_draws dim", print(dim(prior_draws(fn))))
  q("no-prior hyp", print(hypothesis(fn, "x = 0")$hypothesis))
  q("warm default sim", print(fn$fit@sim[c("warmup2", "n_save")]))
  q("update sample_prior", {
    fu <- update(fn, sample_prior = "yes", refresh = 0)
    cat(grep("^prior_", variables(fu), value = TRUE), "\n")
    print(attr(fu$prior, "sample_prior"))
  })
  q("update only", {
    fu <- update(fn, sample_prior = "only", refresh = 0)
    print(posterior_summary(fu, variable = "b_x"))
  })
}
