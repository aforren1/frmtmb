# Item 4: wrong spellings on multivariate fits; extract every
# set_prior() call the message suggests, run it, and report what it
# reaches. Lane arm only.
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2/r2-prelude.R")
source(file.path(R2_ROOT, "dev/mvprior-review2/r2-cases.R"))
d <- r2_data()
ns <- asNamespace("frmtmb")
reach <- function(des, pl) {
  r <- ns$resolve_priorlist(des, pl)
  pt <- des$frame[["par_template"]]
  paste(c(vapply(r$entries, function(e) paste(ns$par_template_names(pt[[e$comp]], e$comp)[e$idx], collapse = "+"), ""),
          if (length(r$lower)) paste0("lb:", paste(names(r$lower), collapse = ",")),
          if (length(r$upper)) paste0("ub:", paste(names(r$upper), collapse = ","))), collapse = " ")
}
extract_calls <- function(msg) {
  m <- gregexpr("set_prior[(]([^()]|[(][^()]*[)])*[)]", msg)
  regmatches(msg, m)[[1]]
}
models <- list(
  mv = list(bf(yg ~ x + (1 | g)) + bf(ym ~ z, sigma ~ x) + set_rescor(FALSE), gaussian()),
  mvnl = list(bf(yg ~ x) + bf(ym ~ a * exp(b * x), a ~ 1 + z, b ~ 1, nl = TRUE) + set_rescor(FALSE), gaussian()),
  mvcat = list(bf(yg ~ x, family = gaussian()) + bf(cat4 ~ x, family = categorical()) + set_rescor(FALSE), NULL),
  mvmix = list(bf(yg ~ x, family = gaussian()) + bf(ym ~ x, family = mixture(gaussian(), gaussian())) + set_rescor(FALSE), NULL)
)
wrong <- list(
  list("mv", set_prior("normal(0, 1)", class = "b")),
  list("mv", set_prior("normal(0, 1)", class = "b", coef = "x")),
  list("mv", set_prior("normal(0, 1)", class = "b", coef = "z")),
  list("mv", set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma")),
  list("mv", set_prior("student_t(3, 0, 2.5)", class = "sigma")),
  list("mv", set_prior("", class = "b", lb = 0)),
  list("mv", set_prior("normal(0, 1)", class = "b", coef = "x", dpar = "sigma")),
  list("mv", set_prior("normal(0, 1)", class = "sd")),
  list("mv", set_prior("normal(0, 1)", class = "b", resp = "Yg")),
  list("mvnl", set_prior("normal(0, 1)", class = "b", nlpar = "a")),
  list("mvnl", set_prior("normal(0, 1)", class = "Intercept", nlpar = "a")),
  list("mvnl", set_prior("normal(0, 1)", class = "b")),
  list("mvnl", set_prior("normal(0, 1)", class = "b", coef = "z", nlpar = "a")),
  list("mvcat", set_prior("normal(0, 1)", class = "b", dpar = "mub")),
  list("mvcat", set_prior("normal(0, 1)", class = "Intercept")),
  list("mvcat", set_prior("normal(0, 1)", class = "b", coef = "x")),
  list("mvmix", set_prior("normal(0, 1)", class = "b")),
  list("mvmix", set_prior("normal(0, 1)", class = "Intercept", dpar = "mu2")),
  list("mvmix", set_prior("student_t(3, 0, 2.5)", class = "sigma1")),
  list("mvmix", set_prior("normal(0, 1)", class = "b", resp = "ym"))
)
des_cache <- list()
for (w in wrong) {
  mn <- w[[1]]; m <- models[[mn]]
  if (is.null(des_cache[[mn]])) {
    des_cache[[mn]] <- tryCatch(ns$prior_design(m[[1]], d, m[[2]], list()), error = function(e) conditionMessage(e))
  }
  des <- des_cache[[mn]]
  s <- w[[2]][[1]]
  f <- c(class = s$class, coef = s$coef, resp = s$resp, dpar = s$dpar, nlpar = s$nlpar)
  cat("\n[", mn, "] ", paste(names(f[nzchar(f)]), f[nzchar(f)], sep = "=", collapse = ", "),
      if (!is.na(s$lb)) " lb=0", "\n", sep = "")
  if (is.character(des)) { cat("  MODEL REFUSED:", des, "\n"); next }
  r <- tryCatch({ reach(des, w[[2]]); "ACCEPTED (not wrong?)" }, error = function(e) conditionMessage(e))
  cat("  msg:", gsub("[[:space:]]+", " ", r), "\n")
  calls <- extract_calls(r)
  if (!length(calls)) { cat("  (no set_prior() call suggested)\n"); next }
  for (cl in calls) {
    rr <- tryCatch(reach(des, eval(parse(text = cl))), error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 150)))
    cat("   ->", cl, "\n      reaches:", rr, "\n")
  }
  both <- tryCatch(reach(des, eval(parse(text = paste(calls, collapse = " + ")))), error = function(e) paste("REFUSED:", substr(conditionMessage(e), 1, 150)))
  cat("   sum of suggestions reaches:", both, "\n")
}
